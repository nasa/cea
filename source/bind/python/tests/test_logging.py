import os
import subprocess
import sys
import textwrap


def _run_python(code: str, *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        check=False,
        env=run_env,
        timeout=20,
    )


def test_default_import_does_not_print_init_paths():
    result = _run_python(
        """
        import cea
        """
    )

    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "Loaded thermo.lib from:" not in combined
    assert "Loaded trans.lib from:" not in combined


def test_info_level_prints_init_paths():
    result = _run_python(
        """
        import cea
        cea.set_log_level(cea.LOG_INFO)
        cea.init()
        """,
        env={"CEA_SKIP_INIT": "1"},
    )

    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "Loaded thermo.lib from:" in result.stdout
    assert "Loaded trans.lib from:" in result.stdout


def test_log_none_suppresses_init_and_backend_logs():
    result = _run_python(
        """
        import cea
        cea.set_log_level(cea.LOG_NONE)
        cea.init()
        _ = cea.Mixture(["H2(L)", "O2(L)"])
        """,
        env={"CEA_SKIP_INIT": "1"},
    )

    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert result.stdout == ""
    assert result.stderr == ""
