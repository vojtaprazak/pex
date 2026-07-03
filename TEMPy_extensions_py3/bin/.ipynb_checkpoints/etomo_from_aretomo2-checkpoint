#!/usr/bin/env python3
"""
make_etomo_directives_from_header.py

Creates per-dataset directives from a template by extracting parameters from the dataset header.

For each dataset directory (e.g. tomog40):
 - detect stack extension (mrc, mrcs, st, ...)
 - call 'header <stack>' (IMOD) and parse pixel spacing and tilt axis angle
 - write directives/<basename>.adoc by replacing the keys:
     setupset.copyarg.name
     setupset.copyarg.pixel
     setupset.copyarg.gold
     setupset.copyarg.rotation
     setupset.copyarg.stackext
     setupset.datasetDirectory
 - optionally run: etomo --headless --directive <adoc>

Requirements:
 - 'header' and 'etomo' must be in PATH (IMOD)
"""

import argparse
import os
import re
import shutil
import subprocess
from glob import glob
import math
from iMOD_comfile import IMOD_comfile

def _parse_linear_xf(path):
    """Return list of (a,b,c,d,tx,ty) for numeric lines in an xf-like file.
       Comments/blank lines are ignored (but could be preserved if needed)."""
    out = []
    with open(path, 'r', encoding='utf-8') as f:
        for lineno, raw in enumerate(f, 1):
            s = raw.strip()
            if not s or s.startswith(('#', ';')):
                continue
            parts = s.split()
            if len(parts) < 6:
                raise ValueError(f"Line {lineno} in {path} has <6 fields: '{s}'")
            try:
                a,b,c,d,tx,ty = map(float, parts[:6])
            except ValueError:
                raise ValueError(f"Line {lineno} in {path} contains non-numeric values: '{s}'")
            out.append((a,b,c,d,tx,ty))
    if not out:
        raise ValueError(f"No numeric transform lines found in {path}")
    return out

def _invert_affine(t):
    a,b,c,d,tx,ty = t
    det = a * d - b * c
    if abs(det) < 1e-12:
        raise ValueError(f"Singular linear part, det={det}: transform {t}")
    ia =  d / det
    ib = -b / det
    ic = -c / det
    id =  a / det
    itx = -(ia * tx + ib * ty)
    ity = -(ic * tx + id * ty)
    return (ia, ib, ic, id, itx, ity)

def _mul_affine(t1, t2):
    """Return t = t1 * t2 (apply t2 then t1)."""
    a1,b1,c1,d1,tx1,ty1 = t1
    a2,b2,c2,d2,tx2,ty2 = t2
    a = a1*a2 + b1*c2
    b = a1*b2 + b1*d2
    c = c1*a2 + d1*c2
    d = c1*b2 + d1*d2
    tx = a1*tx2 + b1*ty2 + tx1
    ty = c1*tx2 + d1*ty2 + ty1
    return (a,b,c,d,tx,ty)

def prexg_to_prexf(prexg_path, out_prexf_path, first_line_identity=True, validate=False, tol=1e-8):
    """
    Convert a linear prexg (g_i) file into a prexf (f_i) file by computing
        f_i = inverse(g_{i-1}) * g_i  (for i>=1)
    and setting f_0 to identity (or to g_0 if first_line_identity=False).

    - prexg_path: path to input prexg (text xf-style file with 6 floats/line)
    - out_prexf_path: path to write prexf
    - first_line_identity: if True, write first f line as 1 0 0 1 0 0; else write f0 = g0
    - validate: if True, recompute g from g0 and f_i and check closeness to input g (useful debug)
    - tol: validation tolerance (max abs diff)
    """
    g = _parse_linear_xf(prexg_path)
    n = len(g)

    f = []
    if first_line_identity:
        f.append((1.0, 0.0, 0.0, 1.0, 0.0, 0.0))
    else:
        f.append(g[0])

    # compute f_i = inv(g_{i-1}) * g_i
    for i in range(1, n):
        inv_prev = _invert_affine(g[i-1])
        fi = _mul_affine(inv_prev, g[i])
        f.append(fi)

    # write file
    with open(out_prexf_path, 'w', encoding='utf-8') as fout:
        for (a,b,c,d,tx,ty) in f:
            fout.write(f"{a:.12g} {b:.12g} {c:.12g} {d:.12g} {tx:.12g} {ty:.12g}\n")

    if validate:
        # reconstruct g_recon by starting from g0 and accumulating f
        g_recon = [g[0]]
        for i in range(1, n):
            g_recon.append(_mul_affine(g_recon[i-1], f[i]))
        # compare
        maxdiff = 0.0
        for i in range(n):
            dif = [abs(g[i][j] - g_recon[i][j]) for j in range(6)]
            maxdiff = max(maxdiff, max(dif))
        if maxdiff > tol:
            raise RuntimeError(f"Validation failed: max diff {maxdiff} > tol {tol}. "
                               "Try first_line_identity=False or inspect transforms.")
    return out_prexf_path

def write_zero_tltxf(xf_path, out_tltxf_path):
     with open(xf_path, 'r', encoding='utf-8') as fin, open(out_tltxf_path, 'w', encoding='utf-8') as fout:
        for line in fin:
            fout.write('1 0 0 1 0 0\n')
            
def xf_to_prexf_shift_only(xf_path, out_prexf_path):
    """
    Read an IMOD xf file (6 floats per non-comment line: a b c d tx ty)
    and write a prexf where we keep only translations (identity 2x2).
    Comment/blank lines are preserved.
    """
    with open(xf_path, 'r', encoding='utf-8') as fin, open(out_prexf_path, 'w', encoding='utf-8') as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith(('#', ';')):
                fout.write(line)
                continue
            parts = s.split()
            if len(parts) < 6:
                # conservative: copy the line if format unexpected
                fout.write(line)
                continue
            try:
                a,b,c,d,tx,ty = parts[:6]
                # write identity linear part, preserve translation
                fout.write(f"1 0 0 1 {tx} {ty}\n")
            except Exception:
                fout.write(line)

# candidate stack extensions to probe (preference order)
KNOWN_STACK_EXTS = ["mrc", "mrcs", "st", "tif", "tiff", "hdf", "h5", "dm4"]

# directive keys to set/replace
KEY_NAME = "setupset.copyarg.name"
KEY_PIXEL = "setupset.copyarg.pixel"
KEY_GOLD = "setupset.copyarg.gold"
KEY_ROTATION = "setupset.copyarg.rotation"
KEY_STACKEXT = "setupset.copyarg.stackext"
KEY_DATASETDIR = "setupset.datasetDirectory"


def run_header_and_capture(stack_path):
    """Run 'header stack_path' and return stdout as text. Raises CalledProcessError if header fails."""
    # 'header' is the IMOD utility that prints stack header information.
    proc = subprocess.run(["header", stack_path], capture_output=True, text=True, check=True)
    return proc.stdout


def parse_header_text(header_text):
    """
    Parse header output for:
      - pixel spacing (first value in 'Pixel spacing (Angstroms)' line)
      - tilt axis angle from 'Tilt axis angle = XX.X' or similar
    Returns (pixel, rotation) where pixel and rotation are strings suitable for writing to directive.
    """
    pixel = None
    rotation = None

    # Pixel spacing line: look for "Pixel spacing (Angstroms)" then numbers
    m_pixel = re.search(r"Pixel spacing.*?([0-9]*\.?[0-9]+)", header_text, re.IGNORECASE)
    if m_pixel:
        pixel = m_pixel.group(1)

    # Tilt axis: many headers show "Tilt axis angle = 84.8, binning = 1 ..." or similar
    m_tilt = re.search(r"Tilt axis angle\s*=\s*([0-9]+\.[0-9]+|[0-9]+)", header_text, re.IGNORECASE)
    if m_tilt:
        rotation = m_tilt.group(1)
    else:
        # alternative pattern sometimes appears after "Titles :" lines
        m_alt = re.search(r"Tilt\s+axis\s+angle\s*=\s*([0-9]+\.[0-9]+|[0-9]+)", header_text, re.IGNORECASE)
        if m_alt:
            rotation = m_alt.group(1)

    return pixel, rotation


def find_stack_ext(dataset_dir, basename, override_ext=None):
    """Return the stack ext (no dot) for the dataset; use override_ext if provided."""
    if override_ext:
        return override_ext.lstrip(".")
    for ext in KNOWN_STACK_EXTS:
        p = os.path.join(dataset_dir, f"{basename}.{ext}")
        if os.path.isfile(p):
            return ext
    # fallback: any file basename.* excluding tlt/xf/adoc
    for path in glob(os.path.join(dataset_dir, f"{basename}.*")):
        if os.path.isfile(path):
            ext = os.path.splitext(path)[1].lstrip(".").lower()
            if ext in ("tlt", "xf", "adoc"):
                continue
            return ext
    return None


def replace_or_append_key(content, key, value):
    """
    Replace an existing line that starts with key (ignoring leading whitespace) with 'key = value'
    or append the line to the end if not present.
    """
    pattern = re.compile(r"(?m)^\s*" + re.escape(key) + r"\s*=.*$")
    new_line = f"{key} = {value}"
    if pattern.search(content):
        content = pattern.sub(new_line, content)
    else:
        if not content.endswith("\n"):
            content += "\n"
        content += "\n" + new_line + "\n"
    return content


def create_per_dataset_adoc(template_path, out_path, replacements):
    with open(template_path, "r", encoding="utf-8") as f:
        tmpl = f.read()
    for key, val in replacements.items():
        tmpl = replace_or_append_key(tmpl, key, val)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(tmpl)


def make_symlink_or_copy(src, dst):
    if os.path.exists(dst):
        return
    try:
        os.symlink(src, dst)
    except (AttributeError, OSError):
        shutil.copy2(src, dst)


def find_and_link_tlt_xf(dataset_dir, basename):
    """Find basename*.tlt and basename*.xf under dataset_dir and link/copy into dataset_dir."""
    for ext in ("tlt", "xf"):
        pattern = os.path.join(dataset_dir, "**", f"{basename}*.{ext}")
        matches = glob(pattern, recursive=True)
        for m in matches:
            dst = os.path.join(dataset_dir, f"{basename}.{ext}")
            if os.path.abspath(m) == os.path.abspath(dst):
                break
            try:
                make_symlink_or_copy(m, dst)
                break
            except Exception:
                continue


def run_etomo_directive(adoc_path, headless=True):
    cmd = ["etomo"]
    if headless:
        cmd.append("--headless")
    cmd += ["--directive", adoc_path]
    print("Running:", " ".join(cmd))
    subprocess.check_call(cmd)


def main():
    p = argparse.ArgumentParser(description="Generate per-dataset IMOD directives from header info and run etomo.")
    p.add_argument("dirs", nargs="+", help="Dataset directories (basename must match tilt-series filename)")
    p.add_argument("--template", "-t", required=True, help="Directive template (.adoc). You can use e.g. the default system templates found in [IMOD directory]/SystemTemplate/cryoSample.adoc ")
    p.add_argument("--gold", type=float, default=10.0, help="Fiducial diameter to set as setupset.copyarg.gold (default: 10)")
    p.add_argument("--stackext", help="Force stack extension (e.g. mrc) rather than auto-detect")
    p.add_argument("--no-run", action="store_true", help="Only create directives; do not run etomo")
    p.add_argument("--skip-existing", action="store_true", help="Skip dataset if directives/<basename>.adoc already exists")
    p.add_argument("--tomo_binning",type=int, default=6)
    p.add_argument("--tomo_thickness",type=int, default=3000)
    p.add_argument("--do_not_run_comfiles", default=False, action='store_true')
    args = p.parse_args()

    cwd = os.getcwd()
    template = os.path.abspath(args.template)
    if not os.path.isfile(template):
        raise SystemExit(f"Template not found: {template}")

    for d in args.dirs:
        if not os.path.isdir(d):
            print(f"Skipping {d}: not a directory")
            continue
        d_abs = os.path.abspath(d)
        basename = os.path.basename(os.path.normpath(d_abs))
        directives_dir = os.path.join(d_abs, "directives")
        os.makedirs(directives_dir, exist_ok=True)
        out_adoc = os.path.join(directives_dir, f"{basename}.adoc")

        if args.skip_existing and os.path.exists(out_adoc):
            print(f"Skipping {d}: {out_adoc} exists")
            continue

        # 1) detect stack extension and path
        stackext = find_stack_ext(d_abs, basename, override_ext=args.stackext)
        if not stackext:
            print(f"ERROR: could not detect a stack file for {basename} in {d_abs}; create {basename}.mrc or specify --stackext")
            continue
        stackfile = os.path.join(d_abs, f"{basename}.{stackext}")
        orig_stack = os.path.join(d_abs, f"{basename}_orig.{stackext}")
        # 2) run header to extract pixel and rotation
        try:
            header_out = run_header_and_capture(stackfile)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: 'header' failed for {stackfile}: {e}; skipping {basename}")
            continue

        pixel, rotation = parse_header_text(header_out)
        if pixel is None:
            print(f"Warning: pixel spacing not found in header for {basename}; leaving blank")
            pixel_val = ""
        else:
            pixel_val = float(pixel)/10

        if rotation is None:
            print(f"Warning: tilt-axis rotation not found in header for {basename}; leaving blank")
            rotation_val = ""
        else:
            rotation_val = rotation

        # 3) ensure tlt/xf are present by linking/copying from subdirs if needed
        find_and_link_tlt_xf(d_abs, basename)

        # 4) assemble replacements (strings)
        replacements = {
            KEY_NAME: basename,
            KEY_PIXEL: pixel_val,
            KEY_GOLD: str(args.gold),
            KEY_ROTATION: rotation_val,
            KEY_STACKEXT: stackext,
            KEY_DATASETDIR: d_abs,
        }

        # 5) create the per-dataset adoc (replace or append keys)
        try:
            create_per_dataset_adoc(template, out_adoc, replacements)
            print(f"Created directive: {out_adoc}")
        except Exception as e:
            print(f"ERROR writing {out_adoc}: {e}")
            continue

        # 6) run etomo with the new directive (unless --no-run)
        if args.no_run:
            continue

        try:
            os.chdir(d_abs)
            run_etomo_directive(out_adoc, headless=True)
            print(f"etomo setup completed for {basename}")
        except subprocess.CalledProcessError as e:
            print(f"etomo failed for {basename}: exit {e.returncode}")
        except FileNotFoundError:
            print("ERROR: 'etomo' not found in PATH. Install IMOD or add etomo to PATH.")
            raise SystemExit(2)
    
        xf_path = os.path.join(d_abs, f"{basename}.xf")
        prexg_path = os.path.join(d_abs, f"{basename}.prexg")
        prexf_path = os.path.join(d_abs, f"{basename}.prexf")
        tltxf_path = os.path.join(d_abs, f"{basename}.tltxf")

        # Preconditions: check xf exists
        if not os.path.isfile(xf_path):
            raise FileNotFoundError(f"{xf_path} not found")

        # create prexf (translation-only)
        shutil.copy2(xf_path, prexg_path)
        write_zero_tltxf(xf_path, tltxf_path)

        prexg_to_prexf(prexg_path, prexf_path)
        if not args.do_not_run_comfiles:
            #eraser
            stack_name = os.path.split(stackfile)[1]
            if os.path.isfile(orig_stack):
                print('%s exists, skipping eraser...' % orig_stack)
            else:
                print('Running eraser...')
                subprocess.check_output('subm eraser.com', shell=True)
                os.rename(stackfile, orig_stack)
                os.rename('%s_fixed.%s' % (basename, stackext), stackfile)
            #prenewst
            subprocess.check_output('subm prenewst.com', shell=True)
            #newst
            print('Generating positioning tomo...')
##            newstcom = IMOD_comfile(d_abs, 'newst.com', make_paths_absolute=False)
##            newstcom.set_val('BinByFactor', args.tomo_binning)
##            if 'SizeToOutputInXandY' in newstcom.dict.keys():
##                del newstcom.dict['SizeToOutputInXandY']
##            newstcom.write_comfile(d_abs)

            newstcom = IMOD_comfile(d_abs, 'newst.com', make_paths_absolute=False)
            newstcom.set_val('BinByFactor', args.tomo_binning)
            # Check if key exists and delete safely using the new API
            if 'SizeToOutputInXandY' in newstcom:
                newstcom.del_val('SizeToOutputInXandY')
            newstcom.write_comfile(d_abs)
            
            subprocess.check_output('subm newst.com', shell=True)
            #tilt
##            tiltcom = IMOD_comfile(d_abs, 'tilt.com', make_paths_absolute=False)
##            tiltcom.dict['IMAGEBINNED'] = args.tomo_binning
##            tiltcom.dict['THICKNESS'] = args.tomo_thickness
##            tiltcom.dict['OutputFile'] = basename + '_full_rec.mrc'
##            del tiltcom.dict['XTILTFILE']
##            tiltcom.write_comfile(d_abs)

            tiltcom = IMOD_comfile(d_abs, 'tilt.com', make_paths_absolute=False)
            tiltcom.set_val('IMAGEBINNED', args.tomo_binning)
            tiltcom.set_val('THICKNESS', args.tomo_thickness)
            tiltcom.set_val('OutputFile', basename + '_full_rec.mrc')

            if 'XTILTFILE' in tiltcom:
                tiltcom.del_val('XTILTFILE')

            tiltcom.write_comfile(d_abs)
            #tiltcom.ScaleShifts
            subprocess.check_output('subm tilt.com', shell=True)

        os.chdir(cwd)

if __name__ == "__main__":
    main()
