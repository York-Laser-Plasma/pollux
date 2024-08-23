import pathlib
import shutil
import subprocess

import pytest


TOP_DIR = pathlib.Path(__file__).parent.parent
GOLDEN_ANSWER_DIR = pathlib.Path(__file__).parent


@pytest.fixture(scope="module")
def build_pollux():
    build_dir = TOP_DIR / "build"
    if not build_dir.exists():
        subprocess.run(f"cmake {TOP_DIR} -B {build_dir}", shell=True, check=True)

    pollux_exe = build_dir / "pollux"
    if not pollux_exe.exists():
        subprocess.run(f"cmake --build {build_dir}", shell=True, check=True)

    return pollux_exe


def test_run(build_pollux, tmp_path):
    shutil.copy(TOP_DIR / "pchris.dat", tmp_path / "pchris.dat")
    subprocess.run(f"{build_pollux}", cwd=tmp_path, check=True)

    original_files = ["orig_fort.11", "orig_fort.12", "orig_fort.13"]
    new_files = ["fort.11", "fort.12", "fort.13"]

    for new, original in zip(new_files, original_files):
        with open(tmp_path / new) as f:
            new_contents = f.read()

        with open(GOLDEN_ANSWER_DIR / original) as f:
            original_contents = f.read()

        assert new_contents == original_contents
