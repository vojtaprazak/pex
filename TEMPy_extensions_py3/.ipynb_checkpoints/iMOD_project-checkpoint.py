"""
# IMOD_project Architecture

This document describes the **architectural rules and invariants** of the `IMOD_project` class after the 2026 refactor. It is written primarily for future maintainers (including future-me) and exists to prevent accidental regression into implicit, stateful, or unsafe behaviour.

This is not user documentation. It explains *why the code is structured the way it is*.

---

## 1. Core design goals

The refactor enforces four non-negotiable goals:

1. **Explicit authority**

   * The project must always know which files are *authoritative*.
   * Authority must never be inferred from filenames or discovery order.

2. **Transactional mutation**

   * Any mutation of IMOD `.com` files must be atomic, reversible, and auditable.
   * Silent overwrites are forbidden.

3. **Separation of intent from execution**

   * Editing `.com` files expresses *intent*.
   * Running IMOD commands produces *effects*.
   * These two steps must never be conflated.

4. **No implicit state coupling**

   * Discovery does not imply correctness.
   * Presence of files does not imply authority.

---

## 2. Project structure

`IMOD_project` is organised into explicit containers. Each container has a single responsibility.

### 2.1 File containers

These containers describe *what exists* on disk. They do not imply authority.

* `image_stacks`

  * `raw`       → `base.mrc`
  * `preali`    → `base_preali.mrc`
  * `ali`       → `base_ali.mrc`

* `volumes`

  * `full_rec`  → `base_full_rec.mrc`
  * `rec`       → `base_rec.mrc` (rotated view)

* `models`

  * fiducials, 3D models, rawfid, etc

* `transforms`

  * `.xf`, `.prexg`, `.prexf`, `local.xf`

* `angles`

  * `.rawtlt`, `.tlt`, `.xtilt`

* `comfiles`

  * `align.com`, `newst.com`, `tilt.com`, etc

* `logfiles`

  * diagnostic output only

These containers are populated by `refresh()`.

---

## 3. Project state (authority)

All *authority* lives in `self.state`.

Key fields:

* `state.current_ali`

  * The authoritative aligned stack.
  * Must equal `image_stacks.ali` when set.

* `state.current_rec`

  * The authoritative reconstruction.
  * Must equal `volumes.full_rec` when set.

* `state.alignment_generation`

* `state.reconstruction_generation`

  * Monotonic counters used for invalidation.

Authority rules:

* A file is authoritative **only if referenced by state**.
* Discovery alone never establishes authority.

---

## 4. Discovery model

### 4.1 `refresh()`

`refresh()` is the **only** method that performs filesystem discovery.

Rules:

* `__init__` must not perform discovery beyond calling `refresh()`.
* `refresh()` populates containers but does not change authority.
* `refresh()` may be called repeatedly and must be idempotent.

---

## 5. Comfile mutation

### 5.1 Single mutation gateway

All `.com` file mutation **must** go through:

```python
with self.modify_comfiles(...) as m:
    ...
```

Direct calls to `IMOD_comfile.write_comfile()` outside this context are forbidden.

### 5.2 Transaction semantics

`modify_comfiles` guarantees:

* Existing `.com` files are backed up to `filename~`.
* Mutations are written atomically.
* On error, backups are restored.
* `refresh()` is called after commit.

---

## 6. Authoritative vs non-authoritative mutation

Not all `.com` mutations imply a change in project truth.

### 6.1 Authoritative mutation (default)

Examples:

* `make_binned_tomo`
* `reconstruct_full`
* (future) `write_ali`

These operations:

* Invalidate downstream authority.
* Increment generation counters.
* Replace canonical outputs.

They use:

```python
with self.modify_comfiles():
```

### 6.2 Non-authoritative mutation

Example:

* `reproject_model`

These operations:

* Temporarily reuse `tilt.com`.
* Must NOT invalidate `current_rec`.
* Must NOT bump generation counters.

They use:

```python
with self.modify_comfiles(authoritative=False):
```

This distinction is intentional and enforced by tests.

---

## 7. Execution model

* Editing `.com` files does **not** execute IMOD.
* Execution is explicit and injected via `_run_comfile`.
* `_run_comfile` is intentionally thin and monkeypatchable for testing.

---

## 8. Invariants

`check_invariants()` enforces:

* Authority consistency (`state` vs containers).
* Generation logic correctness.
* Filesystem sanity for authoritative paths.

It must be called:

* After `refresh()`.
* After any successful `modify_comfiles` commit.
* After authority-creating operations.

---

## 9. What must never happen

The following are architectural violations:

* Writing `.com` files outside `modify_comfiles`.
* Copying `.mrc` or volume files (symlinks only).
* Inferring authority from filenames.
* Mutating `state` based on discovery.
* Using `tilt.com` without explicitly declaring intent.

If any of these occur, the architecture has been broken.

---


"""

import os
import shutil
import warnings
from os.path import join, realpath, isdir, isfile
from subprocess import check_output
import subprocess
from iMOD_comfile import IMOD_comfile
import glob
import shlex
import sys
import tempfile
import numpy as np
from PEETModelParser import PEETmodel
import mrcfile

import os
import shutil
import tempfile
import warnings
from os.path import join, realpath, isfile, exists

# -----------------------------------------------------------------------
# Containers
# -----------------------------------------------------------------------

class ImageStacks(object):
    def __init__(self):
        self.raw = None
        self.preali = None
        self.ali = None

class Volumes(object):
    def __init__(self):
        self.full_rec = None
        self.rec = None

class Models(object):
    def __init__(self):
        self.fid = None
        self.nogaps_fid = None
        self.mod_3d = None
        self.rawfid = None

class Transforms(object):
    def __init__(self):
        self.prexf = None
        self.prexg = None
        self.xf = None
        self.local_xf = None

class Angles(object):
    def __init__(self):
        self.rawtlt = None
        self.tlt = None
        self.xtilt = None

class ComFiles(object):
    def __init__(self):
        self.align = None
        self.prenewst = None
        self.newst = None
        self.tilt = None
        self.ctfcorrection = None

class LogFiles(object):
    def __init__(self):
        self.align = None
        self.newst = None
        self.tilt = None
        self.ctfcorrection = None

class ProjectState(object):
    def __init__(self):
        self.raw_binning = 1
        self.preali_binning = None
        self.ali_binning = None
        self.rec_binning = None
        self.tomo_thickness = None

        self._is_copied = False
        self._is_ctfcorrected = False
        self._is_doseweighted = False

        self.current_ali = None
        self.current_rec = None

        # generation counters for invalidation
        self.alignment_generation = 0
        self.reconstruction_generation = 0

        self.current_ali_confidence = None   # "trusted" | "inherited"
        self.current_rec_confidence = None

# -----------------------------------------------------------------------
# modify_comfiles context manager
# -----------------------------------------------------------------------

class ModifyComfiles(object):
    """
    Context manager to mutate IMOD comfiles safely.

    Usage:
        with ModifyComfiles(proj) as m:
            # m.align is an IMOD_comfile object if present
            m.newst.set_param(block_idx, 'BinByFactor', 6)
            m.tilt.set_param(block_idx, 'OutputFile', 'foo.mrc')

    On exit:
      - original comfiles that would be overwritten are moved to filename~
      - new files are written
      - proj.refresh() is invoked
      - relevant generation counters are incremented
    """

    def __init__(self, proj, authoritative=True):
        self.proj = proj
        self.authoritative = authoritative
        self.root = proj.root
        self._loaded = {}
        self._orig_paths = {}
        self._backups = []
        self._dirty = set()

    def __enter__(self):
        # load known comfiles into memory as IMOD_comfile instances
        from iMOD_comfile import IMOD_comfile

        # map attr name -> filename expected in project root
        mapping = {
            'align': 'align.com',
            'prenewst': 'prenewst.com',
            'newst': 'newst.com',
            'tilt': 'tilt.com',
            'ctfcorrection': 'ctfcorrection.com',
        }

        for attr, fname in mapping.items():
            path = join(self.root, fname)
            self._orig_paths[attr] = path
            if isfile(path):
                # read-only False so that set_param will mark dirty
                try:
                    obj = IMOD_comfile(self.root, fname, make_paths_absolute=False)
                except Exception as e:
                    warnings.warn('[ModifyComfiles] failed to parse %s: %s' % (path, e))
                    obj = None
                self._loaded[attr] = obj
            else:
                self._loaded[attr] = None

        # provide user a simple proxy with attributes
        class Proxy(object):
            pass

        proxy = Proxy()
        for k, v in self._loaded.items():
            setattr(proxy, k, v)
        return proxy

    def _make_backup(self, path):
        if not isfile(path):
            return None
        bak = path + '~'
        # if backup exists, rotate: keep the newest with ~
        try:
            if isfile(bak):
                os.remove(bak)
            os.rename(path, bak)
            self._backups.append((bak, path))
            return bak
        except Exception as e:
            raise RuntimeError('Failed to make backup for %s: %s' % (path, e))

    def _restore_backups(self):
        # try to restore in reverse order
        for bak, orig in reversed(self._backups):
            try:
                if isfile(orig):
                    os.remove(orig)
                os.rename(bak, orig)
            except Exception:
                warnings.warn('[ModifyComfiles] failed to restore backup %s -> %s' % (bak, orig))

    def __exit__(self, exc_type, exc_val, exc_tb):
        # if error inside context, do not write changes, just exit
        if exc_type is not None:
            return False

        # collect comfiles that were modified
        # IMOD_comfile sets _parsed_dirty on mutation via set_param
        to_write = []
        for attr, obj in self._loaded.items():
            if obj is None:
                continue
            dirty = getattr(obj, '_parsed_dirty', False)
            if dirty:
                to_write.append((attr, obj, self._orig_paths.get(attr)))

        if not to_write:
            # nothing to do
            return False

        # perform backups then write; on failure restore backups
        try:
            # 1) backups
            for attr, obj, path in to_write:
                if path and isfile(path):
                    self._make_backup(path)

            # 2) write comfiles
            for attr, obj, path in to_write:
                # write in-place to project root
                try:
                    obj.write_comfile(out_dir=self.root)
                except Exception as e:
                    raise RuntimeError('Failed to write comfile %s: %s' % (path, e))

            # 3) basic validation and state updates
            self._post_commit(to_write)

        except Exception as e:
            # attempt to restore backups and re-raise
            try:
                self._restore_backups()
            finally:
                raise

        return False

    def _post_commit(self, to_write):
        """
        After successful write:
          - refresh project discovery
          - increment generation counters where appropriate
          - sanity-check OutputFile if present
        """
        # simple heuristics: if newst or align modified => alignment_generation++
        # if tilt modified => reconstruction_generation++
        modified_names = [attr for attr, _, _ in to_write]

        # validate presence of OutputFile in written comfiles, warn if missing
        for attr, obj, path in to_write:
            # inspect first executable block for OutputFile param if present
            try:
                for block in getattr(obj, 'blocks', []):
                    params = block.get('params', {})
                    if 'OutputFile' in params:
                        out = params.get('OutputFile')
                        if out is not None:
                            out_path = join(self.root, out)
                            # not necessarily existing yet if execution required
                            if not isfile(out_path):
                                warnings.warn('[ModifyComfiles] OutputFile %s declared in %s but file not present' % (out, path))
                        break
            except Exception:
                warnings.warn('[ModifyComfiles] failed to inspect blocks for %s' % (path))

        # update state counters
        if 'newst' in modified_names or 'align' in modified_names or 'prenewst' in modified_names:
            self.proj.state.alignment_generation += 1
            # if alignment changed, canonical ali/rec become stale
            self.proj.state.current_ali = None
            self.proj.state.current_rec = None

        if self.authoritative and 'tilt' in modified_names:
            self.proj.state.reconstruction_generation += 1
            self.proj.state.current_rec = None

        # refresh discovery
        try:
            self.proj.refresh()
        except Exception as e:
            warnings.warn('[ModifyComfiles] refresh failed after commit: %s' % (e))

        self.proj.check_invariants()

class IMOD_project(object):
    def __init__(self, rec_dir, strict=False):
        self.strict = strict
        self.root = realpath(rec_dir)

        if not os.path.isdir(self.root):
            raise FileNotFoundError(self.root)

        # ------------------------------------------------------------------
        # containers (authoritative structure)
        # ------------------------------------------------------------------
        self.image_stacks = ImageStacks()
        self.volumes = Volumes()
        self.models = Models()
        self.transforms = Transforms()
        self.angles = Angles()
        self.comfiles = ComFiles()
        self.logfiles = LogFiles()
        self.state = ProjectState()
        self.tilt_axis_angle = None
        self.global_excludelist = []

        # ------------------------------------------------------------------
        # identity / provenance
        # ------------------------------------------------------------------
        self.base = None
        self.dataset_name = None
        self.apix = None

        # ------------------------------------------------------------------
        # refresh filesystem discovery
        # ------------------------------------------------------------------
        self.refresh()

    def modify_comfiles(self, authoritative=True):
        return ModifyComfiles(self, authoritative=authoritative)


    def _run_comfile(self, name):
        """
        Execute an IMOD comfile inside the project root.
    
        Uses submfg to preserve IMOD semantics.
        """
        import os
        import subprocess
    
        com_path = os.path.join(self.root, name)
        if not os.path.isfile(com_path):
            raise FileNotFoundError(com_path)
    
        cwd = os.getcwd()
        try:
            os.chdir(self.root)
            subprocess.check_output(
                ["submfg", name],
                stderr=subprocess.STDOUT
            )
        finally:
            os.chdir(cwd)

        
    def _forbid_direct_comfile_write(self):
        raise RuntimeError(
            'Direct comfile writes are forbidden. '
            'Use with self.modify_comfiles(...).'
        )

    def _parse_edf(self, path):
        data = {}
        try:
            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if "=" not in line:
                        continue
                    key, value = line.split("=", 1)
                    data[key.strip()] = value.strip()
        except OSError:
            warnings.warn(f"Failed to read EDF file: {path}")
        return data

    def _warn_if_inherited_authority(self, kind):
        """
        Emit a warning if authoritative input is inherited.
    
        kind: "ali" or "rec"
        """
        import warnings
    
        if kind == "ali":
            conf = self.state.current_ali_confidence
            path = self.state.current_ali
        elif kind == "rec":
            conf = self.state.current_rec_confidence
            path = self.state.current_rec
        else:
            return
    
        if conf == "inherited":
            warnings.warn(
                f"Using inherited authoritative {kind} file: {path}. "
                "Rerun alignment/reconstruction to revalidate."
            )
    
    def _infer_ali_binning(self):
        """
        Infer binning of aligned stack relative to raw.
        Returns int or None.
        """
        # 1) explicit state
        if self.state.ali_binning:
            return self.state.ali_binning
    
        # 2) newst.com
        if self.comfiles.newst and isfile(self.comfiles.newst):
            try:
                com = IMOD_comfile(self.root, 'newst.com', make_paths_absolute=False)
                for block in com.blocks:
                    params = block.get('params', {})
                    for key in ('BinByFactor', 'Binning'):
                        if key in params:
                            return int(params[key])
            except Exception:
                pass
    
        # 3) tilt.com
        if self.comfiles.tilt and isfile(self.comfiles.tilt):
            try:
                com = IMOD_comfile(self.root, 'tilt.com', make_paths_absolute=False)
                for block in com.blocks:
                    params = block.get('params', {})
                    if 'IMAGEBINNED' in params:
                        return int(params['IMAGEBINNED'])
            except Exception:
                pass
    
        return None


    def _infer_preali_binning(self):
        """
        Infer binning of prealigned stack relative to raw.
        Returns int.
        """
        # 1) explicit state
        if self.state.preali_binning:
            return self.state.preali_binning
    
        # 2) prenewst.com
        if self.comfiles.prenewst and isfile(self.comfiles.prenewst):
            try:
                com = IMOD_comfile(self.root, 'prenewst.com', make_paths_absolute=False)
                for block in com.blocks:
                    params = block.get('params', {})
                    for key in ('BinByFactor', 'Binning'):
                        if key in params:
                            return int(params[key])
            except Exception:
                pass
    
        # 3) default: no binning
        return 1
    

    
    def refresh(self):
        root = self.root
    
        # --------------------------------------------------
        # determine base name (non-heuristic, conservative)
        # --------------------------------------------------
        base = None
    
        # 1) EDF metadata (DatasetName, PixelSize, ImageRotationA)
        edf_data = None
        for fname in sorted(os.listdir(root)):
            if fname.lower().endswith('.edf'):
                edf_data = self._parse_edf(join(root, fname))
                if 'Setup.DatasetName' in edf_data:
                    base = edf_data['Setup.DatasetName']
                else:
                    warnings.warn(
                        f"EDF {fname} parsed but Setup.DatasetName not found; "
                        f"keys: {list(edf_data.keys())[:5]}"
                    )
                break

        # --------------------------------------------------
        # tilt axis angle (metadata only)S
        # --------------------------------------------------
        self.tilt_axis_angle = None
        if edf_data is not None:
            val = edf_data.get("Setup.ImageRotationA")
            if val is not None:
                try:
                    self.tilt_axis_angle = float(val)
                except ValueError:
                    warnings.warn(
                        f"Invalid Setup.ImageRotationA value in EDF: {val}"
                    )

    
        # 2) raw stack named base.mrc or base.st (no underscore)
        if base is None:
            raw_candidates = []
            for fname in sorted(os.listdir(root)):
                stem, ext = os.path.splitext(fname)
                if ext.lower() in ('.mrc', '.st') and '_' not in stem:
                    raw_candidates.append(stem)
    
            if len(raw_candidates) == 1:
                base = raw_candidates[0]
    
        # 3) reconstruction fallback
        if base is None:
            for fname in sorted(os.listdir(root)):
                if fname.endswith("_full_rec.mrc"):
                    base = fname[:-len("_full_rec.mrc")]
                    warnings.warn(
                        f"Base name inferred from {fname}; "
                        "no unambiguous dataset identity found"
                    )
                    break
                if fname.endswith("_rec.mrc"):
                    base = fname[:-len("_rec.mrc")]
                    warnings.warn(
                        f"Base name inferred from {fname}; "
                        "no unambiguous dataset identity found"
                    )
                    break
    
            if base is None:
                raise RuntimeError('Unable to determine project base name')
    
        self.base = base

    
        # --------------------------------------------------
        # image stacks (canonical names)
        # --------------------------------------------------
        # self.image_stacks.raw = join(root, '%s.mrc' % base)
        self.image_stacks.raw = None
        self.image_stacks.preali = join(root, '%s_preali.mrc' % base)
        self.image_stacks.ali = join(root, '%s_ali.mrc' % base)
    
        # --------------------------------------------------
        # volumes
        # --------------------------------------------------
        self.volumes.full_rec = join(root, '%s_full_rec.mrc' % base)
        self.volumes.rec = join(root, '%s_rec.mrc' % base)
    
        # --------------------------------------------------
        # models
        # --------------------------------------------------
        self.models.fid = join(root, '%s.fid' % base)
        self.models.nogaps_fid = join(root, '%s_nogaps.fid' % base)
        self.models.mod_3d = join(root, '%s.3dmod' % base)
        self.models.rawfid = join(root, '%s.rawfid' % base)
    
        # --------------------------------------------------
        # transforms
        # --------------------------------------------------
        self.transforms.prexf = join(root, '%s.prexf' % base)
        self.transforms.prexg = join(root, '%s.prexg' % base)
        self.transforms.xf = join(root, '%s.xf' % base)
        self.transforms.local_xf = join(root, '%slocal.xf' % base)
    
        # --------------------------------------------------
        # angles
        # --------------------------------------------------
        self.angles.rawtlt = join(root, '%s.rawtlt' % base)
        self.angles.tlt = join(root, '%s.tlt' % base)
        self.angles.xtilt = join(root, '%s.xtilt' % base)
    
        # --------------------------------------------------
        # comfiles
        # --------------------------------------------------
        self.comfiles.align = join(root, 'align.com')
        self.comfiles.prenewst = join(root, 'prenewst.com')
        self.comfiles.newst = join(root, 'newst.com')
        self.comfiles.tilt = join(root, 'tilt.com')
        self.comfiles.ctfcorrection = join(root, 'ctfcorrection.com')

        self.state.preali_binning = self._infer_preali_binning()
        self.state.ali_binning = self._infer_ali_binning()

        # populate project-level excludelist metadata
        try:
            self.scan_excludelists()
        except Exception:
            # do not fail discovery if scanning fails, just warn
            warnings.warn("scan_excludelists() failed during refresh")
        
        # --------------------------------------------------
        # resolve raw stack from newst.com (supports .mrc, .st, legacy)
        # --------------------------------------------------
        raw = None
        if isfile(self.comfiles.newst):
            try:
                com = IMOD_comfile(self.root, 'newst.com', make_paths_absolute=False)
                for block in com.blocks:
                    params = block.get('params', {})
                    if 'InputFile' in params:
                        raw = join(self.root, params['InputFile'])
                        break
            except Exception as e:
                warnings.warn(f"Failed to parse raw stack from newst.com: {e}")
        
        self.image_stacks.raw = raw

    
        # --------------------------------------------------
        # logfiles
        # --------------------------------------------------
        self.logfiles.align = join(root, 'align.log')
        self.logfiles.newst = join(root, 'newst.log')
        self.logfiles.tilt = join(root, 'tilt.log')
        self.logfiles.ctfcorrection = join(root, 'ctfcorrection.log')

        self.check_invariants()
    
    def check_invariants(self, strict=None):
        """
        Check internal project invariants.
    
        Parameters
        ----------
        strict : bool or None
            If True, raise RuntimeError on violation.
            If False, emit warnings.
            If None, use self.strict.
    
        Returns
        -------
        bool
            True if all invariants hold.
        """
        if strict is None:
            strict = self.strict
    
        errors = []
    
        # --------------------------------------------------
        # authority invariants
        # --------------------------------------------------
        if self.state.current_ali is not None:
            if self.image_stacks.ali is None:
                errors.append('current_ali set but image_stacks.ali is None')
            elif self.state.current_ali != self.image_stacks.ali:
                errors.append(
                    'current_ali mismatch: %s != %s'
                    % (self.state.current_ali, self.image_stacks.ali)
                )
    
        if self.state.current_rec is not None:
            if self.volumes.full_rec is None:
                errors.append('current_rec set but volumes.full_rec is None')
            elif self.state.current_rec != self.volumes.full_rec:
                errors.append(
                    'current_rec mismatch: %s != %s'
                    % (self.state.current_rec, self.volumes.full_rec)
                )
    
        # --------------------------------------------------
        # generation invariants
        # --------------------------------------------------
        if self.state.current_rec is not None:
            if self.state.reconstruction_generation <= 0:
                errors.append(
                    'current_rec set but reconstruction_generation <= 0'
                )
    
        if self.state.current_ali is None and self.state.alignment_generation > 0:
            if self.state.current_rec is not None:
                errors.append(
                    'reconstruction authoritative but alignment invalidated'
                )
    
        # --------------------------------------------------
        # filesystem sanity
        # --------------------------------------------------
        if self.state.current_ali is not None:
            if not os.path.isfile(self.state.current_ali):
                errors.append(
                    'current_ali path does not exist: %s'
                    % self.state.current_ali
                )
    
        if self.state.current_rec is not None:
            if not os.path.isfile(self.state.current_rec):
                errors.append(
                    'current_rec path does not exist: %s'
                    % self.state.current_rec
                )
    
        # --------------------------------------------------
        # report
        # --------------------------------------------------
        if errors:
            msg = 'IMOD_project invariant violation:\n  - ' + '\n  - '.join(errors)
            if strict:
                raise RuntimeError(msg)
            warnings.warn(msg)
            return False
    
        return True

    def copy_project(self, out_root, reset_authority=False, permissive=False):
        """
        Create a new IMOD project rooted in out_root.
    
        Behaviour
        ---------
        - Creates a new project directory named after self.base.
        - Softlinks large files (image_stacks, volumes).
        - Copies small metadata files (comfiles, logs, models, transforms, angles).
        - Instantiates a new IMOD_project.
        - By default, inherits authority and state.
        - reset_authority=True clears authority explicitly.
        - permissive=True allows inherited authority with downgraded confidence.
    
        Returns
        -------
        IMOD_project
            New project instance.
        """
        import os
        import shutil
    
        out_root = os.path.realpath(out_root)
        dest_root = os.path.join(out_root, self.base)
    
        os.makedirs(dest_root, exist_ok=True)
    
        # --------------------------------------------------
        # helper functions
        # --------------------------------------------------
        def safe_symlink(src, dst):
            if not src or not os.path.isfile(src):
                return
            if os.path.exists(dst):
                return
            os.symlink(src, dst)
    
        def safe_copy(src, dst):
            if not src or not os.path.isfile(src):
                return
            if os.path.exists(dst):
                return
            shutil.copy2(src, dst)

        # --------------------------------------------------
        # copy EDF files (metadata, required for identity)
        # --------------------------------------------------
        for fname in os.listdir(self.root):
            if fname.lower().endswith(".edf"):
                src = os.path.join(self.root, fname)
                dst = os.path.join(dest_root, fname)
                safe_copy(src, dst)
                
        # --------------------------------------------------
        # ensure raw stack referenced by newst.com is present
        # --------------------------------------------------
        newst_path = join(self.root, 'newst.com')
        if isfile(newst_path):
            try:
                com = IMOD_comfile(self.root, 'newst.com', make_paths_absolute=False)
                for block in com.blocks:
                    params = block.get('params', {})
                    if 'InputFile' in params:
                        raw = params['InputFile']
                        src = join(self.root, raw)
                        dst = join(dest_root, raw)
                        if isfile(src):
                            safe_symlink(src, dst)
                        break
            except Exception as e:
                warnings.warn(
                    f"Failed to resolve raw stack from newst.com: {e}"
                )
        

        
        # --------------------------------------------------
        # softlink large files
        # --------------------------------------------------
        for path in (
            self.image_stacks.raw,
            self.image_stacks.preali,
            self.image_stacks.ali,
            self.volumes.full_rec,
            self.volumes.rec,
        ):
            if path:
                dst = os.path.join(dest_root, os.path.basename(path))
                safe_symlink(path, dst)
    
        # --------------------------------------------------
        # copy small files
        # --------------------------------------------------
        for container in (
            self.models,
            self.transforms,
            self.angles,
            self.comfiles,
            self.logfiles,
        ):
            for name, path in vars(container).items():
                if not path:
                    continue
                dst = os.path.join(dest_root, os.path.basename(path))
                safe_copy(path, dst)
    
        # --------------------------------------------------
        # instantiate new project
        # --------------------------------------------------
        copied = IMOD_project(dest_root, strict=self.strict)
        copied.refresh()
    
        # --------------------------------------------------
        # transfer or reset state
        # --------------------------------------------------
        if reset_authority:
            copied.state.current_ali = None
            copied.state.current_rec = None
            copied.state.alignment_generation = 0
            copied.state.reconstruction_generation = 0
    
            copied.state.current_ali_confidence = None
            copied.state.current_rec_confidence = None
        else:
            # inherit full state
            # remap authoritative paths to the new project root
            if self.state.current_ali is not None:
                copied.state.current_ali = copied.image_stacks.ali
            else:
                copied.state.current_ali = None
            
            if self.state.current_rec is not None:
                copied.state.current_rec = copied.volumes.full_rec
            else:
                copied.state.current_rec = None

            copied.state.alignment_generation = self.state.alignment_generation
            copied.state.reconstruction_generation = self.state.reconstruction_generation
    
            # confidence handling
            if permissive:
                copied.state.current_ali_confidence = "inherited"
                copied.state.current_rec_confidence = "inherited"
            else:
                copied.state.current_ali_confidence = "trusted"
                copied.state.current_rec_confidence = "trusted"
    
        # --------------------------------------------------
        # final sanity
        # --------------------------------------------------
        if permissive:
            copied.check_invariants(strict=False)
        else:
            copied.check_invariants(strict=True)
    
        return copied
    
    def scan_excludelists(self):
        """
        Scan align.com, newst.com, and tilt.com for excluded views/sections.
    
        Behaviour
        ----------
        - Looks for common parameter names used for exclusion in IMOD comfiles:
          ExcludeList, ExcludeSections, EXCLUDELIST, EXCLUDELIST2
        - Accepts integer ranges and comma separated lists (IMOD style, e.g. "1-8,10,12")
        - Stores sorted unique 1-based indices in self.global_excludelist (as ints)
          so callers can decide how to convert to 0-based if needed.
    
        Returns
        -------
        list[int]
            Sorted unique excluded indices (1-based). Empty list if none found.
        """
        import re
        exclude = set()
    
        def parse_token_list(s):
            """
            Parse strings like "1-4,6,8" and yield ints.
            Accepts ints and ranges, ignores non-parsable tokens.
            """
            s = str(s).strip()
            if not s:
                return []
            out = []
            for tok in re.split(r"[,\s]+", s):
                if not tok:
                    continue
                if "-" in tok:
                    try:
                        a, b = tok.split("-", 1)
                        a = int(a)
                        b = int(b)
                        if b >= a:
                            out.extend(range(a, b + 1))
                    except Exception:
                        continue
                else:
                    try:
                        out.append(int(tok))
                    except Exception:
                        continue
            return out
    
        for comname in ("align.com", "newst.com", "tilt.com"):
            compath = os.path.join(self.root, comname)
            if not os.path.isfile(compath):
                continue
            try:
                com = IMOD_comfile(self.root, comname, make_paths_absolute=False)
            except Exception:
                # fallback: try simple text search for obvious lines
                try:
                    with open(compath, "r") as fh:
                        for line in fh:
                            L = line.strip()
                            if not L:
                                continue
                            # detect Explicit ExcludeList lines
                            if L.lower().startswith("excludelist") or "excludesections" in L.lower():
                                parts = L.split(None, 1)
                                if len(parts) > 1:
                                    for v in parse_token_list(parts[1]):
                                        exclude.add(int(v))
                except Exception:
                    continue
            else:
                # IMOD_comfile parsed OK. Check both com.excludelist if available
                # and the block params for common keys.
                # 1) prefer com.excludelist attribute if present
                if getattr(com, "excludelist", None):
                    for v in com.excludelist:
                        try:
                            exclude.add(int(v))
                        except Exception:
                            pass
    
                # 2) inspect block params for keys containing 'exclude'
                for block in getattr(com, "blocks", []):
                    params = block.get("params", {}) or {}
                    for key, val in params.items():
                        if not isinstance(key, str):
                            continue
                        k = key.lower()
                        if "exclude" not in k:
                            continue
                        # val may be list, string, or None
                        if val is None:
                            continue
                        if isinstance(val, (list, tuple)):
                            for item in val:
                                for v in parse_token_list(item):
                                    exclude.add(int(v))
                        else:
                            for v in parse_token_list(val):
                                exclude.add(int(v))
    
        # Normalise and store
        out = sorted({int(x) for x in exclude if isinstance(x, (int, float))})
        # ensure ints
        out = [int(x) for x in out]
        self.global_excludelist = out
        return out

    
    def make_binned_tomo(self, binning, leave_no_trace=False, allow_regenerate_ali=True):
    
        # must have an ali, or explicitly regenerate it
        if not self.state.current_ali or not isfile(self.state.current_ali):
            if allow_regenerate_ali:
                self.write_ali()
            else:
                raise RuntimeError(
                    "Cannot create binned tomogram: no authoritative aligned stack present"
                )
    
        self._warn_if_inherited_authority("ali")
    
        if not isinstance(binning, int) or binning <= 0:
            raise ValueError("binning must be a positive integer")
    
        base = self.base
        if not base:
            raise RuntimeError("Project base not set; call refresh() first")
    
        # canonical filenames (never change)
        ali_path = join(self.root, '%s_ali.mrc' % base)
        rec_path = join(self.root, '%s_full_rec.mrc' % base)
    
        if not isfile(join(self.root, 'newst.com')) or not isfile(join(self.root, 'tilt.com')):
            raise FileNotFoundError('newst.com and tilt.com must exist in %s' % self.root)
    
        # transactional comfile mutation
        with self.modify_comfiles(authoritative=False) as m:
            newst = m.newst
            tilt = m.tilt
    
            if newst is None or tilt is None:
                raise RuntimeError('newst.com or tilt.com could not be loaded')
    
            # --- newst ---
            newstack_block = None
            for i, block in enumerate(newst.blocks):
                if 'newstack' in block['header'].lower():
                    newstack_block = i
                    break
            if newstack_block is None:
                raise RuntimeError('no newstack block found in newst.com')
    
            params = newst.blocks[newstack_block]['params']
            if 'SizeToOutputInXandY' in params:
                params.pop('SizeToOutputInXandY', None)
                newst.blocks[newstack_block]['separators'].pop('SizeToOutputInXandY', None)
    
            for key in ('BinByFactor', 'Binning'):
                if key in params:
                    newst.set_param(newstack_block, key, binning)
                    break
            else:
                newst.set_param(newstack_block, 'BinByFactor', binning)
    
            # --- tilt ---
            tilt_block = None
            for i, block in enumerate(tilt.blocks):
                if block['header'].lower().startswith('$tilt'):
                    tilt_block = i
                    break
            if tilt_block is None:
                raise RuntimeError('No $tilt block found in tilt.com')
    
            tilt.set_param(tilt_block, 'IMAGEBINNED', binning)
            # DO NOT change InputProjections or OutputFile
    
        # execute
        self._run_comfile('newst.com')
        self._run_comfile('tilt.com')
    
        # authority update (files overwritten in-place)
        if not isfile(ali_path):
            raise RuntimeError('Expected aligned stack not found after newst')
    
        if not isfile(rec_path):
            raise RuntimeError('Expected reconstruction not found after tilt')
    
        self.image_stacks.ali = ali_path
        self.volumes.full_rec = rec_path
    
        self.state.current_ali = ali_path
        self.state.current_rec = rec_path
        self.state.current_ali_confidence = "trusted"
        self.state.current_rec_confidence = "trusted"
        self.state.reconstruction_generation += 1
    
        self.check_invariants(strict=True)
        return rec_path

    
    
    def reconstruct_full(self, regenerate_ali=False, use_shell=False):
        """
        Run reconstruction workflow.
    
        Behaviour:
          - optionally run newst.com to regenerate aligned stack
          - run tilt.com to produce full reconstruction
          - update state.current_rec on success
    
        Parameters
        ----------
        regenerate_ali : bool
            if True, run newst.com before tilt.com
        use_shell : bool
            kept for API compatibility but self._run_comfile uses check_output;
            this flag is ignored here to avoid using shell=True implicitly.
    
        Returns
        -------
        bool
            True if tilt.com executed successfully
        """
        self._warn_if_inherited_authority("ali")
        tilt_path = join(self.root, 'tilt.com')
        if not isfile(tilt_path):
            msg = 'tilt.com not found in %s' % self.root
            if self.strict:
                raise FileNotFoundError(msg)
            warnings.warn(msg)
            return None
    
        if regenerate_ali:
            newst_path = join(self.root, 'newst.com')
            if not isfile(newst_path):
                msg = 'regenerate_ali=True but newst.com not found in %s' % self.root
                if self.strict:
                    raise FileNotFoundError(msg)
                warnings.warn(msg)
            else:
                try:
                    self._run_comfile('newst.com')
                    # after successful newst update state.current_ali if ali exists
                    ali = self.image_stacks.ali
                    if isfile(ali):
                        self.state.current_ali = ali
                        self.state.ali_binning = self.state.ali_binning or None
                except Exception as e:
                    raise RuntimeError('newst.com failed: %s' % e)
    
        # run tilt.com
        try:
            self._run_comfile('tilt.com')
        except Exception as e:
            raise RuntimeError('tilt.com failed: %s' % e)
    
        # set canonical reconstruction if present
        full = self.volumes.full_rec
        if isfile(full):
            self.state.current_rec = full
            self.state.reconstruction_generation += 1
        else:
            warnings.warn('tilt.com finished but %s not found' % full)

        self.state.current_rec_confidence = "trusted"
        self.check_invariants()
        return True
    
    def reproject_model(self, model_3d, out_fid=None, cleanup=False):
        """
        Reproject a 3D model into the aligned stack frame using a
        dedicated reprojection comfile.
        """
    
        if not isfile(model_3d):
            raise FileNotFoundError(model_3d)
    
        if not isfile(self.image_stacks.ali):
            raise FileNotFoundError("Aligned stack not found")
    
        if not isfile(self.comfiles.tilt):
            raise FileNotFoundError("tilt.com not found")
    
        model_3d = os.path.abspath(model_3d)
    
        if out_fid is None:
            base = os.path.splitext(os.path.basename(model_3d))[0]
            out_fid = join(self.root, f"{base}_reprojected.mod")
    
        out_fid = os.path.abspath(out_fid)
    
        tmp_mod = join(self.root, f".tmp_{self.base}_reproj.mod")
        reproj_com = join(self.root, "reproject_model.com")
    
        # --------------------------------------------------
        # 1) transform model into reconstruction frame
        # --------------------------------------------------
        subprocess.check_output(
            ["imodtrans", "-i", self.volumes.full_rec, model_3d, tmp_mod],
            cwd=self.root
        )
    
        # --------------------------------------------------
        # 2) infer binning for reprojection
        # --------------------------------------------------
        binning = self._infer_ali_binning()
        if binning is None:
            raise RuntimeError("Cannot infer IMAGEBINNED for reprojection")
    
        # --------------------------------------------------
        # 3) write reprojection comfile by cloning tilt.com
        # --------------------------------------------------    
        tilt = IMOD_comfile(self.root, 'tilt.com', make_paths_absolute=False)
        
        tilt_block = None
        for i, block in enumerate(tilt.blocks):
            if block['header'].lower().startswith('$tilt'):
                tilt_block = i
                break
        
        if tilt_block is None:
            raise RuntimeError("No $tilt block found in tilt.com")
        
        # override ONLY what is necessary
        tilt.set_param(tilt_block, 'ProjectModel', os.path.basename(tmp_mod))
        tilt.set_param(tilt_block, 'OutputFile', os.path.basename(out_fid))
        
        # write as a separate comfile
        tilt.write_comfile(
            out_dir=self.root,
            change_name='reproject_model.com'
        )

    
        # --------------------------------------------------
        # 4) run reprojection
        # --------------------------------------------------
        self._run_comfile("reproject_model.com")
    
        if not isfile(out_fid):
            raise RuntimeError("Reprojection failed")
    
        # --------------------------------------------------
        # cleanup
        # --------------------------------------------------
        if cleanup:
            for p in (tmp_mod, reproj_com):
                if isfile(p):
                    os.remove(p)
    
        self.check_invariants()
        return out_fid
    
            
    def write_ali(self, ali_output_path=None, binning=None,
                  use_shell=False, honour_excludelist=False):
        """
        Produce the aligned stack by running the newst.com newstack block.
    
        Behaviour (migrated)
        --------------------
        - Uses the modify_comfiles() transactional API to mutate newst.com safely.
        - Runs newst.com (via _run_comfile) to generate the aligned file.
        - On success updates self.state.current_ali and increments alignment_generation.
        - Does not copy or move any .mrc files (no heavy file operations).
        - Caller can supply ali_output_path (absolute or relative to project root).
        - Caller can supply binning (positive int) to override binning in newst.com.
        - honour_excludelist: if True, inject current global_excludelist into newst.com.
    
        Returns
        -------
        str
            Absolute path to produced aligned stack (inside project root).
        """
        import os
    
        # sanity checks
        if not os.path.isfile(os.path.join(self.root, "newst.com")):
            raise FileNotFoundError("newst.com")
    
        raw_stack = self.image_stacks.raw
        if not raw_stack or not os.path.isfile(raw_stack):
            raise FileNotFoundError(raw_stack or "raw_stack not set")
    
        # We'll mutate newst.com via the transactional API
        with self.modify_comfiles() as m:
            # m.newst is an IMOD_comfile instance for newst.com (make_paths_absolute=False)
            newst = getattr(m, "newst", None)
            if newst is None:
                raise RuntimeError("modify_comfiles() did not provide newst.com")
    
            # find the newstack block index (raise if not present)
            newstack_block = None
            for i, block in enumerate(newst.blocks):
                if "newstack" in block["header"].lower():
                    newstack_block = i
                    break
            if newstack_block is None:
                raise RuntimeError("no newstack block found in newst.com")
    
            params = newst.blocks[newstack_block]["params"]
    
            if "InputFile" not in params:
                raise RuntimeError("newst.com missing InputFile")
    
            # honour_excludelist -> set ExcludeSections if available
            if honour_excludelist and getattr(self, "global_excludelist", None):
                newst.set_param(
                    newstack_block,
                    "ExcludeSections",
                    self.global_excludelist,
                    separator=","
                )
    
            input_stack = params["InputFile"]
            # ensure the input stack exists relative to project root
            if not os.path.isfile(os.path.join(self.root, input_stack)):
                raise FileNotFoundError(input_stack)
    
            # check any referenced transform files exist (defensive)
            for key in ("TransformFile", "XfFile", "XFFile"):
                if key in params and not os.path.isfile(os.path.join(self.root, params[key])):
                    raise FileNotFoundError(params[key])
    
            # choose output path and set OutputFile in comfile
            if ali_output_path is not None:
                out_path_rel = ali_output_path
                # if absolute, write relative name into comfile but ultimately run from self.root
                if os.path.isabs(out_path_rel):
                    # write absolute path into comfile only if caller supplied absolute
                    newst.set_param(newstack_block, "OutputFile", out_path_rel)
                    out_path_abs = out_path_rel
                else:
                    newst.set_param(newstack_block, "OutputFile", out_path_rel)
                    out_path_abs = os.path.join(self.root, out_path_rel)
            elif "OutputFile" in params:
                out_path_rel = params["OutputFile"]
                out_path_abs = os.path.join(self.root, out_path_rel)
            else:
                base, _ = os.path.splitext(os.path.basename(input_stack))
                out_path_rel = base + "_ali.mrc"
                newst.set_param(newstack_block, "OutputFile", out_path_rel)
                out_path_abs = os.path.join(self.root, out_path_rel)
    
            # handle binning override
            if binning is not None:
                if not isinstance(binning, int) or binning <= 0:
                    raise ValueError("binning must be a positive integer")
    
                # remove explicit output size if present (it would conflict)
                if "SizeToOutputInXandY" in params:
                    # mutate in-memory representation to avoid writing contradictory fields
                    params.pop("SizeToOutputInXandY", None)
                    newst.blocks[newstack_block]["separators"].pop("SizeToOutputInXandY", None)
    
                # set BinByFactor/Binning parameter
                for key in ("BinByFactor", "Binning"):
                    if key in params:
                        newst.set_param(newstack_block, key, binning)
                        break
                else:
                    newst.set_param(newstack_block, "BinByFactor", binning)
    
            # At this point, the comfile has been mutated in-memory by the context manager.
            # The context manager will write the comfile when exiting the 'with' block.
            # We need to exit the with-block to let it write, then execute newst.com.
            # But we also want to ensure the run happens with the newly-written comfile.
            # So we store the out_path_abs and ask the context manager to commit.
            # Exiting the `with` will have written newst.com into self.root.
    
        # Outside the with-block: comfile has been written (or modified transactionally).
        # Execute the command (this can be monkeypatched in tests via _run_comfile)
        try:
            self._run_comfile("newst.com")
        except subprocess.CalledProcessError as e:
            # _run_comfile may raise; report clearly
            raise RuntimeError(f"newst.com execution failed: {e}")
    
        # sanity check: aligned output exists
        if not os.path.isfile(out_path_abs):
            # If newst.com wrote a different path (e.g. absolute), try resolving from comfile
            # reload discovery minimally
            self.refresh()
            if not os.path.isfile(out_path_abs):
                raise RuntimeError(f"newst.com completed but output not found: {out_path_abs}")
    
        # update authoritative state now that a new aligned stack exists
        self.state.current_ali = out_path_abs
        # bump alignment generation so downstream state invalidates recs etc
        self.state.alignment_generation = (getattr(self.state, "alignment_generation", 0) or 0) + 1
        self.state.current_ali_confidence = "trusted"
        self.check_invariants()
        
        return out_path_abs
        
    def backtransform_model_to_raw(self, model_path, out_model=None, binning=None):
        """
        Backtransform a model from aligned/reconstruction coordinates
        into raw stack coordinates.
    
        Parameters
        ----------
        model_path : str
            Path to input model (picked on ali or rec)
        out_model : str or None
            Output model path. Default: <model>_raw.mod
        binning : int or None
            Binning factor of the aligned stack relative to raw.
            If None, inferred from state.
    
        Returns
        -------
        str
            Absolute path to raw-space model
        """
    
        if not isfile(model_path):
            raise FileNotFoundError(model_path)
    
        if not self.transforms.xf or not isfile(self.transforms.xf):
            raise FileNotFoundError("Alignment transform (.xf) not found")
    
        if not self.image_stacks.raw or not isfile(self.image_stacks.raw):
            raise FileNotFoundError("Raw stack not found")
    
        model_path = os.path.abspath(model_path)
    
        if out_model is None:
            base, _ = os.path.splitext(model_path)
            out_model = base + "_raw.mod"
    
        out_model = os.path.abspath(out_model)
    
        # infer binning (permissive)
        if binning is None:
            binning = self._infer_ali_binning()
            if binning is None:
                raise RuntimeError(
                    "Failed to infer binning from state or comfiles"
                )
    
        # IMOD expects scale = 1/binning
        scale = 1.0 / float(binning)
    
        cmd = [
            "xfmodel",
            "-back",
            "-scale", str(scale),
            "-xform",
            self.transforms.xf,
            model_path,
            out_model,
        ]
        #print(' '.join(cmd))
        subprocess.check_output(cmd, cwd=self.root)
    
        if not isfile(out_model):
            raise RuntimeError("Backtransform failed, output not created")
    
        return out_model







