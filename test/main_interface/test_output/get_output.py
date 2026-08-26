import os
import shutil
import subprocess
from pathlib import Path


def resolve_cea_exe(repo_root: Path) -> str:
    cea_exe = os.environ.get("CEA_EXE")
    if cea_exe is not None:
        return cea_exe

    candidate_exes = []
    executable_names = ("cea.exe", "cea") if os.name == "nt" else ("cea",)

    run_dir = os.environ.get("CEA_RUN_DIR")
    if run_dir is not None:
        candidate_exes.extend(
            Path(run_dir).expanduser() / name for name in executable_names
        )

    # Prefer binaries built from this checkout to avoid stale external builds.
    for name in executable_names:
        candidate_exes.extend(sorted(repo_root.glob(f"build*/source/{name}")))
        candidate_exes.extend(sorted(repo_root.glob(f"build*/**/source/{name}")))

    # Keep the historical default location as a final fallback.
    candidate_exes.extend(
        Path(f"~/git/cea/build-dev/source/{name}").expanduser()
        for name in executable_names
    )

    existing_exes = [path for path in candidate_exes if path.exists()]
    if existing_exes:
        return str(max(existing_exes, key=lambda path: path.stat().st_mtime))

    cea_exe = shutil.which("cea")
    if cea_exe is not None:
        return cea_exe

    raise FileNotFoundError(
        "Could not locate the `cea` executable. Set CEA_EXE/CEA_RUN_DIR or build CEA in this repository."
    )


def run_tests(test_names):

    script_dir = Path(__file__).resolve().parent.parent
    repo_root = script_dir.parent.parent
    cea_exe = resolve_cea_exe(repo_root)

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
