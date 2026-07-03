import os
from iMOD_comfile import IMOD_comfile


HERE = os.path.dirname(__file__)
FILES_DIR = os.path.join(HERE, "files")


def _read(path):
    with open(path, "r") as f:
        return f.read()


def _write(path, txt):
    with open(path, "w") as f:
        f.write(txt)


def test_all_comfiles_roundtrip_identity(tmp_path):
    """
    For every .com file in tests/iMOD_comfile/files:
      - copy to tmp_path
      - parse with IMOD_comfile
      - write back without modification
      - assert bytewise identical output

    This is the strongest possible invariant.
    """
    files = [f for f in os.listdir(FILES_DIR) if f.endswith(".com")]
    assert len(files) > 0, "No .com files found for testing"

    for i in range(len(files)):
        fname = files[i]
        src = os.path.join(FILES_DIR, fname)
        dst = os.path.join(tmp_path, fname)

        original = _read(src)
        _write(dst, original)

        obj = IMOD_comfile(tmp_path, fname, make_paths_absolute=False)
        obj.write_comfile(tmp_path, change_name=fname)

        rewritten = _read(dst)

        assert original == rewritten, "Roundtrip mismatch for %s" % fname


def test_all_comfiles_mutate_and_rewrite(tmp_path):
    """
    For every .com file:
      - parse it
      - apply conservative, generic mutations
      - write to new file
      - ensure output is non-empty and contains valid structure

    This test does NOT require semantic correctness of parameters,
    only that the library can safely rewrite reasonable changes.
    """
    files = [f for f in os.listdir(FILES_DIR) if f.endswith(".com")]
    assert len(files) > 0, "No .com files found for testing"

    for i in range(len(files)):
        fname = files[i]
        src = os.path.join(FILES_DIR, fname)
        dst = os.path.join(tmp_path, fname)

        _write(dst, _read(src))

        obj = IMOD_comfile(tmp_path, fname, make_paths_absolute=False)

        # Apply generic mutations to first block only
        if obj.blocks:
            blk = obj.blocks[0]

            # 1) If RADIAL exists, tweak it
            if "RADIAL" in blk["params"]:
                obj.set_param(0, "RADIAL", 0.2, separator="")

            # 2) If EXCLUDELIST exists, normalize it
            for k in ["EXCLUDELIST", "ExcludeList", "ExcludeSections"]:
                if k in blk["params"]:
                    obj.set_param(0, k, [1, 3, 5], separator=",")

            # 3) Otherwise, add a harmless numeric parameter
            if "TestParameter" not in blk["params"]:
                obj.set_param(0, "TestParameter", 1)

        out_name = "mutated_" + fname
        obj.write_comfile(tmp_path, change_name=out_name)

        out_path = os.path.join(tmp_path, out_name)
        assert os.path.isfile(out_path), "Failed to write %s" % out_name

        out_txt = _read(out_path)

        # Minimal sanity checks
        assert len(out_txt.strip()) > 0
        assert "$" in out_txt, "No command header found in %s" % out_name
