import pytest
import numpy as np
from PEETModelParser import PEETmodel


def test_empty_model_creates_without_error():
    m = PEETmodel()
    assert len(m) == 0


def test_add_point_increments_length():
    m = PEETmodel()
    m.add_point(0, 0, [10.0, 20.0, 30.0])
    assert len(m) == 1
    m.add_point(0, 0, [11.0, 21.0, 31.0])
    assert len(m) == 2


def test_read_two_points_mod(data_dir):
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    assert m.no_of_obj >= 1
    assert len(m) >= 2


def test_get_vector_returns_vector(data_dir):
    from Vector import Vector
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    v = m.get_vector(0)
    assert isinstance(v, Vector)


def test_distance_is_positive(data_dir):
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    d = m.distance(0, 1)
    assert d > 0.0


def test_round_trip(data_dir, tmp_path):
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    original_points = m.get_all_points().copy()

    out = str(tmp_path / 'roundtrip.mod')
    m.write_model(out)

    m2 = PEETmodel(out)
    assert len(m2) == len(m)
    np.testing.assert_allclose(m2.get_all_points(), original_points, atol=1e-3)
