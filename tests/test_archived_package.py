import importlib


def test_legacy_kartik_eli_is_archived():
    try:
        importlib.import_module('pergamon.archive.kartik_eli')
    except ImportError:
        return
    raise AssertionError('pergamon.archive.kartik_eli should be excluded from the supported API and fail on import.')
