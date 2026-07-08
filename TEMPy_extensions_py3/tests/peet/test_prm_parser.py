import pytest
from PEETPRMParser import PEETPRMFile


def test_read_run1_prm_fnoutput(data_dir):
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    assert prm.prm_dict['fnOutput'] == 'r'


def test_read_run1_prm_fnmodparticle_is_set(data_dir):
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    assert prm.prm_dict.get('fnModParticle') is not None


def test_round_trip(data_dir, tmp_path):
    src = str(data_dir / 'run1' / 'r.prm')
    prm = PEETPRMFile(src)

    out = str(tmp_path / 'roundtrip.prm')
    prm.write_prm_file(out)

    prm2 = PEETPRMFile(out)
    assert prm2.prm_dict['fnOutput'] == prm.prm_dict['fnOutput']
    assert prm2.prm_dict.get('fnModParticle') == prm.prm_dict.get('fnModParticle')
