import os
from pathlib import Path

from pergamon.paths import retr_pathgaia
from pergamon.vetting_workflow import retr_pathvetting


def test_retr_pathvetting_creates_expected_directories(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path / 'pergamon'))

    dictpath = retr_pathvetting(strgextn='SyntheticPopulation_TESS')

    assert os.path.isdir(dictpath['pathbase'])
    assert os.path.isdir(dictpath['pathdatapipe'])
    assert os.path.isdir(dictpath['pathvisupipe'])
    assert os.path.isdir(dictpath['pathpopl'])
    assert os.path.isdir(dictpath['pathvisucnfg'])
    assert os.path.isdir(dictpath['pathdatacnfg'])
    assert Path(dictpath['pathbasepergamon']) == tmp_path / 'pergamon'
    assert os.path.normpath(dictpath['pathpopl']).endswith(os.path.normpath('SyntheticPopulation_TESS'))


def test_retr_pathgaia_creates_expected_directories(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path / 'pergamon'))

    dictpath = retr_pathgaia()

    assert os.path.isdir(dictpath['pathbase'])
    assert os.path.isdir(dictpath['pathdata'])
    assert os.path.isdir(dictpath['pathimag'])
    assert os.path.normpath(dictpath['pathbase']).endswith(os.path.normpath('pergamon'))