import os
import numpy as np
from iMOD_comfile import IMOD_comfile


def test_tilt_parse_roundtrip_and_command(tmp_path):
    here = os.path.dirname(__file__)
    files_dir = os.path.join(here, "files")

    src = os.path.join(files_dir, "tilt.com")
    dst = os.path.join(tmp_path, "tilt.com")

    with open(src, "r") as f:
        content = f.read()

    with open(dst, "w") as f:
        f.write(content)

    obj = IMOD_comfile(tmp_path, "tilt.com", make_paths_absolute=False)

    # Parse
    assert obj.dict["InputFile"] == "test.st"
    assert obj.dict["OutputFile"] == "test.rec"

    excl = np.atleast_1d(obj.dict["EXCLUDELIST"])
    assert set(excl) == {1, 3, 4, 5}

    # Round trip
    obj.out_dir = tmp_path
    obj.write_comfile()

    with open(dst, "r") as f:
        written = f.read()

    assert "InputFile" in written
    assert "OutputFile" in written
    assert "RADIAL" in written
    assert "EXCLUDELIST" in written

    # Command conversion
    cmd = obj.get_command_list()

    assert cmd[0] == "tilt"
    assert "-InputFile" in cmd
    assert "-OutputFile" in cmd

    # Default excluded parameter
    assert "-RADIAL" not in cmd

    # Excludelist passed
    assert "-EXCLUDELIST" in cmd
    excl_val = cmd[cmd.index("-EXCLUDELIST") + 1]
    for x in ["1", "3", "4", "5"]:
        assert x in excl_val
