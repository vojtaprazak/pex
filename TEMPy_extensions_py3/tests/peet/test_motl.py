import pytest
import numpy as np
from PEETMotiveList import PEETMotiveList


def test_empty_motl_creates_without_error():
    m = PEETMotiveList()
    assert len(m) == 0


def test_add_particle_increments_length():
    m = PEETMotiveList()
    m.add_empty_pcle(angles=[10.0, 20.0, 30.0])
    assert len(m) == 1
    m.add_empty_pcle(angles=[0.0, 0.0, 0.0])
    assert len(m) == 2


def test_read_iter1_particle_count(data_dir):
    m = PEETMotiveList(str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv'))
    assert len(m) == 258


def test_read_iter2_particle_count(data_dir):
    m = PEETMotiveList(str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter2.csv'))
    assert len(m) == 258


def test_get_angles_returns_three_floats(data_dir):
    m = PEETMotiveList(str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv'))
    angles = m.get_angles_by_list_index(0)
    assert len(angles) == 3
    assert all(isinstance(a, float) for a in angles)


def test_round_trip(data_dir, tmp_path):
    src = str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv')
    m = PEETMotiveList(src)
    out = str(tmp_path / 'roundtrip.csv')
    m.write_PEET_motive_list(out)

    m2 = PEETMotiveList(out)
    assert len(m2) == len(m)
    for i in range(len(m)):
        a1 = m.get_angles_by_list_index(i)
        a2 = m2.get_angles_by_list_index(i)
        np.testing.assert_allclose(a1, a2, atol=1e-3)


def test_angles_to_norm_vec_are_unit_vectors(data_dir):
    m = PEETMotiveList(str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv'))
    nv = m.angles_to_norm_vec()
    for v in nv:
        length = (v.x**2 + v.y**2 + v.z**2) ** 0.5
        assert abs(length - 1.0) < 1e-6, f"norm_vec has length {length}, expected 1.0"
