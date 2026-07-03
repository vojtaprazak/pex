import os
import shutil
import pytest
from os.path import join, isfile
import subprocess 

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------

def write_file(path, text):
    with open(path, 'w') as f:
        f.write(text)

def read_file(path):
    with open(path, 'r') as f:
        return f.read()

def minimal_project(tmp_path):
    """
    Create a minimal fake IMOD project directory with:
      - base.mrc
      - base.rawtlt
      - newst.com
      - tilt.com
    No IMOD binaries required.
    """
    root = tmp_path / 'proj'
    root.mkdir()

    base = 'test'

    write_file(root / ('%s.mrc' % base), 'RAWSTACK')
    write_file(root / ('%s.rawtlt' % base), '0\n1\n2\n')

    write_file(
        root / 'newst.com',
        '$newstack\nInputFile %s.mrc\nOutputFile %s_ali.mrc\n' % (base, base)
    )

    write_file(
        root / 'tilt.com',
        '$tilt\nInputProjections %s_ali.mrc\nOutputFile %s_full_rec.mrc\n' % (base, base)
    )

    return root, base


# ---------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------

def test_refresh_populates_containers(tmp_path):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    # image stacks
    assert proj.image_stacks.raw.endswith('%s.mrc' % base)
    assert proj.image_stacks.preali.endswith('%s_preali.mrc' % base)
    assert proj.image_stacks.ali.endswith('%s_ali.mrc' % base)

    # volumes
    assert proj.volumes.full_rec.endswith('%s_full_rec.mrc' % base)
    assert proj.volumes.rec.endswith('%s_rec.mrc' % base)

    # containers exist
    assert proj.models is not None
    assert proj.transforms is not None
    assert proj.angles is not None
    assert proj.comfiles is not None
    assert proj.logfiles is not None
    assert proj.state is not None


def test_modify_comfiles_creates_backup(tmp_path):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    orig = read_file(join(root, 'newst.com'))

    with proj.modify_comfiles() as m:
        newst = m.newst
        assert newst is not None

        # find first block
        idx = 0
        newst.set_param(idx, 'OutputFile', '%s_b2_ali.mrc' % base)

    assert isfile(join(root, 'newst.com~'))
    assert isfile(join(root, 'newst.com'))

    new = read_file(join(root, 'newst.com'))
    bak = read_file(join(root, 'newst.com~'))

    assert bak == orig
    assert new != orig


def test_modify_comfiles_rollback_on_exception(tmp_path):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    orig = read_file(join(root, 'newst.com'))

    with pytest.raises(RuntimeError):
        with proj.modify_comfiles() as m:
            newst = m.newst
            newst.set_param(0, 'OutputFile', 'BROKEN.mrc')
            raise RuntimeError('boom')

    assert read_file(join(root, 'newst.com')) == orig
    assert not isfile(join(root, 'newst.com~'))


def test_alignment_invalidates_reconstruction_state(tmp_path):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    proj.state.current_ali = join(root, '%s_ali.mrc' % base)
    proj.state.current_rec = join(root, '%s_full_rec.mrc' % base)
    proj.state.alignment_generation = 0

    with proj.modify_comfiles() as m:
        m.newst.set_param(0, 'OutputFile', '%s_b4_ali.mrc' % base)

    assert proj.state.current_ali is None
    assert proj.state.current_rec is None
    assert proj.state.alignment_generation == 1


def test_reconstruct_sets_current_rec(tmp_path, monkeypatch):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    # fake execution: touching output file
    def fake_run(name):
        if name == 'tilt.com':
            write_file(join(root, '%s_full_rec.mrc' % base), 'REC')

    monkeypatch.setattr(proj, '_run_comfile', fake_run)

    ok = proj.reconstruct_full()
    assert ok is True
    assert proj.state.current_rec.endswith('%s_full_rec.mrc' % base)
    assert proj.state.reconstruction_generation == 1


def test_no_mrc_copying(tmp_path, monkeypatch):
    from iMOD_project import IMOD_project

    root, base = minimal_project(tmp_path)
    proj = IMOD_project(str(root))
    proj.refresh()

    # fail test if any .mrc file is copied
    def forbid_copy(src, dst, *a, **k):
        if src.endswith('.mrc'):
            raise AssertionError('Attempted to copy large MRC file')
        return shutil.copy2(src, dst)

    monkeypatch.setattr(shutil, 'copy2', forbid_copy)

    # only mutate comfiles, do not execute
    with proj.modify_comfiles() as m:
        m.newst.set_param(0, 'OutputFile', '%s_b8_ali.mrc' % base)


def test_reproject_model_does_not_change_authority(tmp_path, monkeypatch):
    from iMOD_project import IMOD_project
    import subprocess
    import mrcfile

    root, base = minimal_project(tmp_path)

    # fake reconstruction so method passes checks
    write_file(root / ('%s_full_rec.mrc' % base), 'REC')

    # fake input model
    write_file(root / 'dummy.mod', 'MOD')

    proj = IMOD_project(str(root))
    proj.refresh()

    proj.state.current_rec = proj.volumes.full_rec
    proj.state.reconstruction_generation = 1

    # --------------------------------------------------
    # stub external commands
    # --------------------------------------------------
    def fake_check_output(cmd, *a, **k):
        # simulate model2point writing model.xyz
        if isinstance(cmd, (list, tuple)) and 'model2point' in cmd[0]:
            # last argument is xyz_path
            xyz_path = cmd[-1]
            with open(xyz_path, 'w') as f:
                f.write('10 10 0\n')
        return None
    
    monkeypatch.setattr(
        subprocess,
        'check_output',
        fake_check_output
    )


    monkeypatch.setattr(
        proj,
        '_run_comfile',
        lambda name: write_file(
            root / 'test_reprojected.3dfid',
            'FID'
        )
    )

    # --------------------------------------------------
    # stub mrcfile.open to avoid real MRC parsing
    # --------------------------------------------------
    class FakeMrc(object):
        def __init__(self):
            class H(object):
                ny = 100
            self.header = H()
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        mrcfile,
        'open',
        lambda *a, **k: FakeMrc()
    )

    out = proj.reproject_model(
        str(root / 'dummy.mod'),
        out_fid=str(root / 'test_reprojected.3dfid')
    )

    assert out.endswith('test_reprojected.3dfid')
    assert proj.state.current_rec == proj.volumes.full_rec

def test_write_ali_sets_current_ali_and_generation(tmp_path, monkeypatch):
    from iMOD_project import IMOD_project
    import shutil

    # --------------------------------------------------
    # minimal project
    # --------------------------------------------------
    root, base = minimal_project(tmp_path)

    # raw stack must exist for write_ali
    # minimal_project already created base.mrc

    proj = IMOD_project(str(root))
    proj.refresh()

    # initial state
    assert proj.state.current_ali is None
    assert proj.state.alignment_generation == 0

    # --------------------------------------------------
    # stub execution: pretend newst.com ran and produced ali
    # --------------------------------------------------
    ali_path = root / ('%s_ali.mrc' % base)

    def fake_run_comfile(name):
        # simulate newst.com creating aligned stack
        if name == "newst.com":
            with open(ali_path, "w") as f:
                f.write("ALI")

    monkeypatch.setattr(proj, "_run_comfile", fake_run_comfile)

    # --------------------------------------------------
    # forbid copying of .mrc files
    # --------------------------------------------------
    def forbid_copy(src, dst, *a, **k):
        if src.endswith(".mrc"):
            raise AssertionError("write_ali attempted to copy an MRC file")
        return shutil.copy2(src, dst)

    monkeypatch.setattr(shutil, "copy2", forbid_copy)

    # --------------------------------------------------
    # run write_ali
    # --------------------------------------------------
    out = proj.write_ali()

    # --------------------------------------------------
    # assertions
    # --------------------------------------------------
    assert out.endswith('%s_ali.mrc' % base)
    assert proj.state.current_ali == str(ali_path)
    assert proj.state.alignment_generation == 1

def test_copy_project_permissive_inherits_authority_with_confidence(tmp_path):
    from iMOD_project import IMOD_project

    # --------------------------------------------------
    # original project
    # --------------------------------------------------
    root, base = minimal_project(tmp_path)

    # simulate authoritative alignment + reconstruction
    write_file(root / f"{base}_ali.mrc", "ALI")
    write_file(root / f"{base}_full_rec.mrc", "REC")

    proj = IMOD_project(str(root))
    proj.refresh()

    proj.state.current_ali = proj.image_stacks.ali
    proj.state.current_rec = proj.volumes.full_rec
    proj.state.alignment_generation = 1
    proj.state.reconstruction_generation = 1
    proj.state.current_ali_confidence = "trusted"
    proj.state.current_rec_confidence = "trusted"

    proj.check_invariants(strict=True)

    # --------------------------------------------------
    # permissive copy
    # --------------------------------------------------
    out_root = tmp_path / "copies"
    copied = proj.copy_project(
        str(out_root),
        reset_authority=False,
        permissive=True
    )

    # --------------------------------------------------
    # assertions
    # --------------------------------------------------
    assert copied.state.current_ali == copied.image_stacks.ali
    assert copied.state.current_rec == copied.volumes.full_rec

    assert copied.state.current_ali_confidence == "inherited"
    assert copied.state.current_rec_confidence == "inherited"

    # permissive mode should not raise
    copied.check_invariants(strict=False)

