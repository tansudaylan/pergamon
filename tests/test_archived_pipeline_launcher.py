import runpy
from pathlib import Path


def test_north_pipeline_delegates_with_northern_sectors(monkeypatch, capsys):
    launcher_path = (
        Path(__file__).parents[1] / "pergamon" / "archive" / "kartik_eli" / "mf_pipeline_north.py"
    )
    original_run_path = runpy.run_path
    namespace = original_run_path(str(launcher_path), run_name="mf_pipeline_north_test")
    calls = []
    monkeypatch.setattr(runpy, "run_path", lambda *args, **kwargs: calls.append((args, kwargs)))

    namespace["main"]()

    pipeline_path = launcher_path.with_name("mf_pipeline.py")
    assert calls == [((str(pipeline_path),), {"init_globals": {"SECTORS": range(15, 29)}})]
    assert capsys.readouterr().out == f"Reading from {pipeline_path}...\n"