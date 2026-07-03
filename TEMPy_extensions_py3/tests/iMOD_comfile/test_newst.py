import os
from iMOD_comfile import IMOD_comfile


def test_newst_parse_roundtrip_and_command(tmp_path):
    here = os.path.dirname(__file__)
    files_dir = os.path.join(here, "files")

    src = os.path.join(files_dir, "newst.com")
    dst = os.path.join(tmp_path, "newst.com")

    with open(src, "r") as f:
        content = f.read()

    with open(dst, "w") as f:
        f.write(content)

    obj = IMOD_comfile(tmp_path, "newst.com", make_paths_absolute=False)

    # Parse
    assert obj.dict["InputFile"] == "test.mrc"
    assert obj.dict["OutputFile"] == "test.st"
    assert obj.dict["TransformFile"] == "test.xf"

    # Round trip
    obj.out_dir = tmp_path
    obj.write_comfile()

    with open(dst, "r") as f:
        written = f.read()

    assert "InputFile" in written
    assert "OutputFile" in written
    assert "TransformFile" in written

    # Command conversion
    cmd = obj.get_command_list()

    assert cmd[0] == "newstack"
    assert "-InputFile" in cmd
    assert cmd[cmd.index("-InputFile") + 1] == "test.mrc"
    assert "-OutputFile" in cmd
    assert cmd[cmd.index("-OutputFile") + 1] == "test.st"
    assert "-TransformFile" in cmd
    assert cmd[cmd.index("-TransformFile") + 1] == "test.xf"
