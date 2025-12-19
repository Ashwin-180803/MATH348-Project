import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_mean_curvature


def assert_mean_curvature_equal(H, expected, tol=1e-12):
    """
    Compare computed mean curvature H with expected value.

    Handles:
      - Exact symbolic comparisons
      - Numeric expressions with a tolerance

    Parameters:
        H : sympy expression
        expected : sympy expression or numeric
        tol : float, tolerance for numeric comparison
    """
    H = sp.sympify(H)
    expected = sp.sympify(expected)

    diff = sp.simplify(H - expected)
    if diff == 0:
        return True

    free_syms = list(H.free_symbols.union(expected.free_symbols))
    subs = {s: 1 for s in free_syms} if free_syms else {}
    H_val = float(H.subs(subs).evalf())
    expected_val = float(expected.subs(subs).evalf())

    if abs(H_val - expected_val) > tol:
        raise AssertionError(
            f"Mean curvature mismatch:\nComputed: {H}\nExpected: {expected}\nNumeric diff: {H_val - expected_val}"
        )


# ======================================================================
# 1. PLANE
# H = 0
# ======================================================================


def test_mean_curvature_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 2. CYLINDER
# H = -1/(2a)
# ======================================================================


def test_mean_curvature_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -1 / (2 * R)
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -1 / (2 * a)
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 3. SPHERE
# H = -1/a
# ======================================================================


def test_mean_curvature_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -1 / R
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -1 / a
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 4. TORUS
# H = -cos(v)/(2(a+b cos v)) - 1/(2b)
# ======================================================================


def test_mean_curvature_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -sp.cos(v) / (2 * (R + r * sp.cos(v))) - 1 / (2 * r)
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = -sp.cos(v) / (2 * (a + b * sp.cos(v))) - 1 / (2 * b)
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 5. PARABOLOID
# H = (1 + 4 u^2 + 4 v^2)^(-3/2) * (1 + 2(u^2+v^2)) simplified
# ======================================================================


def test_mean_curvature_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 2 / ((1 + 4 * u**2 + 4 * v**2) ** (3 / 2)) * (1 + 2 * u**2 + 2 * v**2)
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = (a * (1 + 4 * b**2 * v**2) + b * (1 + 4 * a**2 * u**2)) / (
        1 + 4 * a**2 * u**2 + 4 * b**2 * v**2
    ) ** (3 / 2)

    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 6. HYPERBOLIC PARABOLOID
# H = 0
# ======================================================================


def test_mean_curvature_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 7. ELLIPTIC TYPE SURFACE
# H computed from trace(I^-1 II)/2
# ======================================================================


def test_mean_curvature_elliptic_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u + v, u - v, u * v]
    Xu = sp.Matrix([1, 1, v])
    Xv = sp.Matrix([1, -1, u])
    n = (Xu.cross(Xv)).normalized()
    Xuu = sp.Matrix([0, 0, 0])
    Xuv = sp.Matrix([0, 0, 1])
    Xvv = sp.Matrix([0, 0, 0])
    L = Xuu.dot(n)
    M = Xuv.dot(n)
    N = Xvv.dot(n)
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).trace() / 2

    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    n = (Xu.cross(Xv)).normalized()
    Xuu = sp.Matrix([0, 0, 0])
    Xuv = sp.Matrix([0, 0, a * b])
    Xvv = sp.Matrix([0, 0, 0])
    L = Xuu.dot(n)
    M = Xuv.dot(n)
    N = Xvv.dot(n)
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).trace() / 2

    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 8. CATENOID
# H = 0
# ======================================================================


def test_mean_curvature_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 9. HELICOID
# H = 0
# ======================================================================


def test_mean_curvature_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    expected = 0
    assert_mean_curvature_equal(H, expected)


# ======================================================================
# 10. NONLINEAR SURFACE
# H computed from trace(I^-1 II)/2
# ======================================================================


def test_mean_curvature_weird_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    Xu = sp.Matrix([2 * u * v, v**2, sp.exp(u + v)])
    Xv = sp.Matrix([u**2, 2 * u * v, sp.exp(u + v)])
    Xuu = sp.Matrix([2 * v, 0, sp.exp(u + v)])
    Xuv = sp.Matrix([2 * u, 2 * v, sp.exp(u + v)])
    Xvv = sp.Matrix([0, 2 * u, sp.exp(u + v)])
    n = (Xu.cross(Xv)).normalized()
    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).trace() / 2

    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    assert_mean_curvature_equal(H, expected)


def test_mean_curvature_weird_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    Xu = sp.Matrix([2 * a * u * v, b * v**2, a * sp.exp(a * u + b * v)])
    Xv = sp.Matrix([a * u**2, 2 * b * u * v, b * sp.exp(a * u + b * v)])
    Xuu = sp.Matrix([2 * a * v, 0, a**2 * sp.exp(a * u + b * v)])
    Xuv = sp.Matrix([2 * a * u, 2 * b * v, a * b * sp.exp(a * u + b * v)])
    Xvv = sp.Matrix([0, 2 * b * u, b**2 * sp.exp(a * u + b * v)])
    n = (Xu.cross(Xv)).normalized()
    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).trace() / 2

    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    H, _ = compute_mean_curvature(X, [u, v])
    H, _ = compute_mean_curvature(X, [u, v])
    assert_mean_curvature_equal(H, expected)
