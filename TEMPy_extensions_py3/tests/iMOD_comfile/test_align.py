import os
from iMOD_comfile import IMOD_comfile


def test_align_parse_roundtrip_and_command(tmp_path):
    here = os.path.dirname(__file__)
    files_dir = os.path.join(here, "files")

    src = os.path.join(files_dir, "align.com")
    dst = os.path.join(tmp_path, "align.com")

    with open(src, "r") as f:
        content = f.read()

    with open(dst, "w") as f:
        f.write(content)

    obj = IMOD_comfile(tmp_path, "align.com", make_paths_absolute=False)

    # Parse first block
    assert obj.dict["InputFile"] == "tilt.st"
    assert obj.dict["OutputFile"] == "align.xf"

    # Parse second block
    assert obj.b_dict["InputFile"] == "align.xf"
    assert obj.b_dict["OutputFile"] == "final.xf"

    # Round trip
    obj.out_dir = tmp_path
    obj.write_comfile()

    with open(dst, "r") as f:
        written = f.read()

    assert "$tiltalign -StandardInput" in written
    assert "$xfproduct -StandardInput" in written

    # Command conversion uses first program only
    cmd = obj.get_command_list()

    assert cmd[0] == "tiltalign"
    assert "-InputFile" in cmd
    assert cmd[cmd.index("-InputFile") + 1] == "tilt.st"
    assert "final.xf" not in cmd
