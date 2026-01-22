import os
import shutil
import subprocess
from pathlib import Path

def run_tests(test_names):

    script_dir = Path(__file__).resolve().parent.parent
    repo_root = script_dir.parent.parent
    run_dir = Path(os.environ.get("CEA_RUN_DIR", "~/git/cea/build-dev/source")).expanduser()
    default_exe = run_dir / "cea"
    cea_exe = os.environ.get("CEA_EXE")
    if cea_exe is None:
        if default_exe.exists():
            cea_exe = str(default_exe)
        else:
            cea_exe = shutil.which("cea")
    if cea_exe is None:
        raise FileNotFoundError("Could not locate `cea` executable. Set CEA_EXE or CEA_RUN_DIR, or add `cea` to PATH.")

    for test in test_names:
        # Execute the code on the input file
        print(f"Running {test}")
        input_base = script_dir / test
        subprocess.run([cea_exe, str(input_base)], cwd=repo_root, check=True)
        print()

    return

if __name__=="__main__":
    # test_names = ["example1", "example2", "example3",
    #               "example4","example5", "example6",
    #               "example7", "example8","example9",
    #               "example10","example11", "example12",
    #               "example13", "example14"]
    test_names = ["example1"]
    run_tests(test_names)
