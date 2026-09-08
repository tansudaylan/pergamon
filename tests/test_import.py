import importlib

import pergamon


def test_import_pergamon_package():
    assert hasattr(pergamon, '__file__')


def test_main_module_exposes_init():
    main = importlib.import_module('pergamon.main')
    assert hasattr(main, 'init')
