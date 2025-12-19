import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_christoffel_symbols


def expected_christoffel(E, F, G, u, v):
    E_u, E_v = sp.diff(E, u), sp.diff(E, v)
    F_u, F_v = sp.diff(F, u), sp.diff(F, v)
    G_u, G_v = sp.diff(G, u), sp.diff(G, v)
    I = sp.Matrix([[E, F], [F, G]])
    I_inv = I.inv()
    gamma_sub_uu = sp.Matrix([[0.5 * E_u], [F_u - 0.5 * E_v]])
    gamma_sub_uv = sp.Matrix([[0.5 * E_v], [0.5 * G_u]])
    gamma_sub_vv = sp.Matrix([[F_v - 0.5 * G_u], [0.5 * G_v]])
    return {
        "gamma_u_sub_uu": (I_inv.row(0) * gamma_sub_uu)[0],
        "gamma_v_sub_uu": (I_inv.row(1) * gamma_sub_uu)[0],
        "gamma_u_sub_uv": (I_inv.row(0) * gamma_sub_uv)[0],
        "gamma_v_sub_uv": (I_inv.row(1) * gamma_sub_uv)[0],
        "gamma_u_sub_vv": (I_inv.row(0) * gamma_sub_vv)[0],
        "gamma_v_sub_vv": (I_inv.row(1) * gamma_sub_vv)[0],
    }


def assert_christoffel_symbols_equal(gamma, expected, tol=1e-12):
    """
    Compare computed Christoffel symbols with expected symbols.

    Parameters:
        gamma : dict of sympy expressions
        expected : dict of sympy expressions
        tol : float, numeric tolerance
    """
    for key in expected:
        computed = sp.sympify(gamma[key])
        exp = sp.sympify(expected[key])
        diff = sp.simplify(computed - exp)
        if diff == 0:
            continue

        free_syms = list(computed.free_symbols.union(exp.free_symbols))
        subs = {s: 1 for s in free_syms} if free_syms else {}
        c_val = float(computed.subs(subs).evalf())
        e_val = float(exp.subs(subs).evalf())
        if abs(c_val - e_val) > tol:
            raise AssertionError(
                f"Christoffel symbol {key} mismatch:\nComputed: {computed}\nExpected: {exp}\nNumeric diff: {c_val - e_val}"
            )


# ======================================================================
# 1. PLANE
# ======================================================================
def test_christoffel_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    expected = {key: 0 for key in gamma}
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    expected = {key: 0 for key in gamma}
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 2. CYLINDER
# X(u,v) = (R cos u, R sin u, b v)
# ======================================================================
def test_christoffel_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    # Compute E,F,G
    E = R**2
    F = 0
    G = 1
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    E = a**2
    F = 0
    G = b**2
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 3. SPHERE
# X(u,v) = (a sin u cos v, a sin u sin v, a cos u)
# ======================================================================
def test_christoffel_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    E = R**2
    F = 0
    G = R**2 * sp.sin(u) ** 2
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    E = a**2
    F = 0
    G = a**2 * sp.sin(u) ** 2
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 4. TORUS
# X(u,v) = ((a + b cos v)cos u, (a + b cos v) sin u, b sin v)
# ======================================================================
def test_christoffel_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix(
        [-(R + r * sp.cos(v)) * sp.sin(u), (R + r * sp.cos(v)) * sp.cos(u), 0]
    )
    Xv = sp.Matrix(
        [-r * sp.sin(v) * sp.cos(u), -r * sp.sin(v) * sp.sin(u), r * sp.cos(v)]
    )
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix(
        [-(a + b * sp.cos(v)) * sp.sin(u), (a + b * sp.cos(v)) * sp.cos(u), 0]
    )
    Xv = sp.Matrix(
        [-b * sp.sin(v) * sp.cos(u), -b * sp.sin(v) * sp.sin(u), b * sp.cos(v)]
    )
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 5. PARABOLOID
# X(u,v) = (u,v,u^2+v^2)
# ======================================================================
def test_christoffel_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, 0, 2 * u])
    Xv = sp.Matrix([0, 1, 2 * v])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, 0, 2 * a * u])
    Xv = sp.Matrix([0, 1, 2 * b * v])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 6. HYPERBOLIC PARABOLOID
# X(u,v) = (u,v,u^2-v^2)
# ======================================================================
def test_christoffel_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, 0, 2 * u])
    Xv = sp.Matrix([0, 1, -2 * v])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, 0, 2 * a * u])
    Xv = sp.Matrix([0, 1, -2 * b * v])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 7. ELLIPTIC SURFACE
# X(u,v) = (u+v, a u-b v, a b u v)
# ======================================================================
def test_christoffel_elliptic_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u + v, u - v, u * v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, 1, v])
    Xv = sp.Matrix([1, -1, u])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 8. CATENOID
# X(u,v) = (a cosh v cos u, a cosh v sin u, b v)
# ======================================================================
def test_christoffel_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([-sp.cosh(v) * sp.sin(u), sp.cosh(v) * sp.cos(u), 0])
    Xv = sp.Matrix([sp.sinh(v) * sp.cos(u), sp.sinh(v) * sp.sin(u), 1])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([-a * sp.cosh(v) * sp.sin(u), a * sp.cosh(v) * sp.cos(u), 0])
    Xv = sp.Matrix([a * sp.sinh(v) * sp.cos(u), a * sp.sinh(v) * sp.sin(u), b])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 9. HELICOID
# X(u,v) = (a u cos v, a u sin v, b v)
# ======================================================================
def test_christoffel_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([sp.cos(v), sp.sin(v), 0])
    Xv = sp.Matrix([-u * sp.sin(v), u * sp.cos(v), 1])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([a * sp.cos(v), a * sp.sin(v), 0])
    Xv = sp.Matrix([-a * u * sp.sin(v), a * u * sp.cos(v), b])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


# ======================================================================
# 10. NONLINEAR SURFACE
# X(u,v) = (a u^2 v, b u v^2, exp(a u + b v))
# ======================================================================
def test_christoffel_weird_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([2 * u * v, v**2, sp.exp(u + v)])
    Xv = sp.Matrix([u**2, 2 * u * v, sp.exp(u + v)])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)


def test_christoffel_weird_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True)
    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    gamma, _ = compute_christoffel_symbols(X, [u, v])
    Xu = sp.Matrix([2 * a * u * v, b * v**2, a * sp.exp(a * u + b * v)])
    Xv = sp.Matrix([a * u**2, 2 * b * u * v, b * sp.exp(a * u + b * v)])
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    expected = expected_christoffel(E, F, G, u, v)
    assert_christoffel_symbols_equal(gamma, expected)
