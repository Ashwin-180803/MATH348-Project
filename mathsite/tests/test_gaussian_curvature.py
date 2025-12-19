import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_gaussian_curvature
import sympy as sp


def assert_gaussian_curvature_equal(K, expected, tol=1e-12):
    """
    Compare computed Gaussian curvature K with expected value.

    Handles:
      - Exact symbolic comparisons
      - Numeric expressions with tolerance

    Parameters:
        K : sympy expression
        expected : sympy expression or numeric
        tol : float
    """
    K = sp.sympify(K)
    expected = sp.sympify(expected)

    diff = sp.simplify(K - expected)
    if diff == 0:
        return True

    free_syms = list(K.free_symbols.union(expected.free_symbols))
    subs = {s: 1 for s in free_syms} if free_syms else {}
    K_val = float(K.subs(subs).evalf())
    expected_val = float(expected.subs(subs).evalf())

    if abs(K_val - expected_val) > tol:
        raise AssertionError(
            f"Gaussian curvature mismatch:\nComputed: {K}\nExpected: {expected}\nNumeric diff: {K_val - expected_val}"
        )


# ======================================================================
# 1. PLANE
# K = 0
# ======================================================================


def test_gaussian_curvature_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 0
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 0
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 2. CYLINDER
# K = 0 (principal curvatures 1/a, 0)
# ======================================================================


def test_gaussian_curvature_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 0
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 0
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 3. SPHERE
# K = 1/a^2
# ======================================================================


def test_gaussian_curvature_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 1 / R**2
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 1 / a**2
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 4. TORUS
# K = (1/b)*(cos(v)/(a+b cos(v))) (product of principal curvatures)
# ======================================================================


def test_gaussian_curvature_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 1 / r * (sp.cos(v) / (R + r * sp.cos(v)))
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 1 / b * (sp.cos(v) / (a + b * sp.cos(v)))
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 5. PARABOLOID
# K = 4 / (1 + 4u^2 + 4v^2)^2
# ======================================================================


def test_gaussian_curvature_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 4 / (1 + 4 * u**2 + 4 * v**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = 4 * a * b / (1 + 4 * a**2 * u**2 + 4 * b**2 * v**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 6. HYPERBOLIC PARABOLOID
# K = -4 / (1 + 4u^2 + 4v^2)^2
# ======================================================================


def test_gaussian_curvature_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -4 / (1 + 4 * u**2 + 4 * v**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -4 * a * b / (1 + 4 * a**2 * u**2 + 4 * b**2 * v**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 7. ELLIPTIC SURFACE
# K = det(I^-1 II)
# ======================================================================


def test_gaussian_curvature_elliptic_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u + v, u - v, u * v]
    Xu = sp.Matrix([1, 1, v])
    Xv = sp.Matrix([1, -1, u])
    n = (Xu.cross(Xv)).normalized()
    Xuu = sp.Matrix([0, 0, 0])
    Xuv = sp.Matrix([0, 0, 1])
    Xvv = sp.Matrix([0, 0, 0])
    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).det()

    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])

    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    n = (Xu.cross(Xv)).normalized()
    Xuu = sp.Matrix([0, 0, 0])
    Xuv = sp.Matrix([0, 0, a * b])
    Xvv = sp.Matrix([0, 0, 0])
    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = (I.inv() * II).det()

    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])

    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 8. CATENOID
# K = -1/(a^2 cosh^4 v)
# ======================================================================


def test_gaussian_curvature_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -1 / sp.cosh(v) ** 4
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -(a**2) * b**2 / (a**2 * sp.sinh(v) ** 2 + b**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 9. HELICOID
# K = -b^2 / (a^2 u^2 + b^2)^2
# ======================================================================


def test_gaussian_curvature_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -1 / (1 + u**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    expected = -(b**2) / (a**2 * u**2 + b**2) ** 2
    assert_gaussian_curvature_equal(K, expected)


# ======================================================================
# 10. NONLINEAR SURFACE
# K = det(I^-1 II)
# ======================================================================


def test_gaussian_curvature_weird_surface_numeric():
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
    expected = (I.inv() * II).det()

    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    assert_gaussian_curvature_equal(K, expected)


def test_gaussian_curvature_weird_surface_symbolic():
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
    expected = (I.inv() * II).det()

    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    K, _ = compute_gaussian_curvature(X, [u, v])
    K, _ = compute_gaussian_curvature(X, [u, v])
    assert_gaussian_curvature_equal(K, expected)
