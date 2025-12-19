import sympy as sp
import pytest
from calcapp.utils.functions import compute_frenet_serret_apparatus


def assert_vectors_equal(v1, v2):
    assert len(v1) == len(v2)
    for a, b in zip(v1, v2):
        assert sp.simplify(a - b) == 0


def assert_scalars_equal(a, b):
    assert sp.simplify(a - b) == 0


def assert_nan_or_zero(x):
    """Accept scalar that is either exactly 0 or NaN."""
    if x.is_nan:
        return
    assert sp.simplify(x) == 0


def assert_vector_nan_or_zero(v):
    """Accept vector whose components are each 0 or NaN."""
    for comp in v:
        if comp.is_nan:
            continue
        assert sp.simplify(comp) == 0


# ===============================================================
# 1. LINE: X(t) = (t, 2t, 0)
# ===============================================================


def test_frenet_line_numeric():
    t = sp.symbols("t", real=True)
    X = [t, 2 * t, 0]

    with pytest.raises(
        ValueError, match="Frenet-Serret apparatus undefined for straight lines."
    ):
        compute_frenet_serret_apparatus(X, t)


def test_frenet_line_symbolic():
    t = sp.symbols("t", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * t + b, c * t + d, 0]

    with pytest.raises(
        ValueError, match="Frenet-Serret apparatus undefined for straight lines."
    ):
        compute_frenet_serret_apparatus(X, t)


# ===============================================================
# 2. CIRCLE: X(t) = (R cos t, R sin t, 0)
# ===============================================================


def test_frenet_circle_numeric():
    t = sp.symbols("t", real=True)
    R = 2.7

    X = [R * sp.cos(t), R * sp.sin(t), 0]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    T_expected = sp.Matrix([-sp.sin(t), sp.cos(t), 0])
    N_expected = sp.Matrix([-sp.cos(t), -sp.sin(t), 0])
    B_expected = sp.Matrix([0, 0, 1])
    kappa = 1 / R

    assert_vectors_equal(res["T"], T_expected)
    assert_vectors_equal(res["N"], N_expected)
    assert_vectors_equal(res["B"], B_expected)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], 0)


def test_frenet_circle_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)

    X = [a * sp.cos(t), a * sp.sin(t), 0]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    T_expected = sp.Matrix([-sp.sin(t), sp.cos(t), 0])
    N_expected = sp.Matrix([-sp.cos(t), -sp.sin(t), 0])
    B_expected = sp.Matrix([0, 0, 1])
    kappa = 1 / a

    assert_vectors_equal(res["T"], T_expected)
    assert_vectors_equal(res["N"], N_expected)
    assert_vectors_equal(res["B"], B_expected)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], 0)


# ===============================================================
# 3. SPHERE MERIDIAN: X(t) = (R sin t, 0, R cos t)
# ===============================================================


def test_frenet_sphere_meridian_numeric():
    t = sp.symbols("t", real=True)
    R = 1.8

    X = [R * sp.sin(t), 0, R * sp.cos(t)]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    T_expected = sp.Matrix([sp.cos(t), 0, -sp.sin(t)])
    N_expected = sp.Matrix([-sp.sin(t), 0, -sp.cos(t)])
    B_expected = sp.Matrix([0, 1, 0])
    kappa = 1 / R

    assert_vectors_equal(res["T"], T_expected)
    assert_vectors_equal(res["N"], N_expected)
    assert_vectors_equal(res["B"], B_expected)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], 0)


def test_frenet_sphere_meridian_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)

    X = [a * sp.sin(t), 0, a * sp.cos(t)]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    T_expected = sp.Matrix([sp.cos(t), 0, -sp.sin(t)])
    N_expected = sp.Matrix([-sp.sin(t), 0, -sp.cos(t)])
    B_expected = sp.Matrix([0, 1, 0])
    kappa = 1 / a

    assert_vectors_equal(res["T"], T_expected)
    assert_vectors_equal(res["N"], N_expected)
    assert_vectors_equal(res["B"], B_expected)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], 0)


# ===============================================================
# 4. HELIX: X(t) = (cos t, sin t, b t)
# ===============================================================


def test_frenet_helix_numeric():
    t = sp.symbols("t", real=True)
    b = 2

    X = [sp.cos(t), sp.sin(t), b * t]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    v = sp.sqrt(1 + b**2)

    T = sp.Matrix([-sp.sin(t) / v, sp.cos(t) / v, b / v])
    N = sp.Matrix([-sp.cos(t), -sp.sin(t), 0])
    B = T.cross(N)
    kappa = 1 / (1 + b**2)
    tau = b / (1 + b**2)

    assert_vectors_equal(res["T"], T)
    assert_vectors_equal(res["N"], N)
    assert_vectors_equal(res["B"], B)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], tau)


def test_frenet_helix_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)

    X = [a * sp.cos(t), a * sp.sin(t), b * t]
    res, _ = compute_frenet_serret_apparatus(X, t)
    res, _ = compute_frenet_serret_apparatus(X, t)

    v = sp.sqrt(a**2 + b**2)

    T = sp.Matrix([-a * sp.sin(t) / v, a * sp.cos(t) / v, b / v])
    N = sp.Matrix([-sp.cos(t), -sp.sin(t), 0])
    B = T.cross(N)
    kappa = a / (a**2 + b**2)
    tau = b / (a**2 + b**2)

    assert_vectors_equal(res["T"], T)
    assert_vectors_equal(res["N"], N)
    assert_vectors_equal(res["B"], B)
    assert_scalars_equal(res["kappa"], kappa)
    assert_scalars_equal(res["tau"], tau)
