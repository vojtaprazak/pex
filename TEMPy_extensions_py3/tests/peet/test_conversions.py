import pytest
import numpy as np
from conversions import euler_from_matrix, zyz_to_zxz, zxz_to_zyz
from transformations import euler_matrix


# --- regression: our compat wrapper must accept 3x3 ---

def test_euler_from_matrix_accepts_3x3():
    """Regression: PyPI transformations rejects 3x3; our wrapper must pad it."""
    mat3 = np.eye(3)
    result = euler_from_matrix(mat3, 'rzxz')
    assert len(result) == 3


def test_euler_from_matrix_accepts_4x4():
    mat4 = np.eye(4)
    result = euler_from_matrix(mat4, 'rzxz')
    assert len(result) == 3


def test_euler_from_matrix_3x3_rotation_gives_sensible_angles():
    """A 90-degree rotation around Z should give a known ZXZ decomposition."""
    from transformations import rotation_matrix
    mat4 = rotation_matrix(np.radians(90), [0, 0, 1])
    mat3 = mat4[:3, :3]
    result = euler_from_matrix(mat3, 'rzxz')
    assert len(result) == 3
    assert all(np.isfinite(a) for a in result)


# --- Euler angle conversion round-trips ---

@pytest.mark.parametrize("z1,x,z2", [
    (0.0,  0.0,   0.0),
    (30.0, 45.0,  60.0),
    (90.0, 10.0, -45.0),
    (180.0, 90.0, 270.0),
])
def test_zxz_to_zyz_roundtrip(z1, x, z2):
    zyz = zxz_to_zyz(z1, x, z2)
    back = zyz_to_zxz(*zyz)
    # Reconstruct rotation matrices and compare — angle decompositions
    # are not unique, so compare via matrix product instead.
    m1 = euler_matrix(np.radians(z1), np.radians(x), np.radians(z2), 'rzxz')[:3, :3]
    m2 = euler_matrix(np.radians(back[0]), np.radians(back[1]), np.radians(back[2]), 'rzxz')[:3, :3]
    np.testing.assert_allclose(m1, m2, atol=1e-10)


@pytest.mark.parametrize("z1,y,z2", [
    (0.0,  0.0,  0.0),
    (30.0, 45.0, 60.0),
    (90.0, 10.0, -45.0),
])
def test_zyz_to_zxz_roundtrip(z1, y, z2):
    zxz = zyz_to_zxz(z1, y, z2)
    back = zxz_to_zyz(*zxz)
    m1 = euler_matrix(np.radians(z1), np.radians(y), np.radians(z2), 'rzyz')[:3, :3]
    m2 = euler_matrix(np.radians(back[0]), np.radians(back[1]), np.radians(back[2]), 'rzyz')[:3, :3]
    np.testing.assert_allclose(m1, m2, atol=1e-10)
