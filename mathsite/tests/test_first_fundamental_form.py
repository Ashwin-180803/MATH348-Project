import sympy as sp
from calcapp.utils.parser import parse_input
import pytest
from calcapp.utils.functions import compute_first_fundamental_form


def assert_matrix_equal(A, B):
    assert A.shape == B.shape

    symsA = sorted(A.free_symbols, key=lambda s: s.name)
    symsB = sorted(B.free_symbols, key=lambda s: s.name)

    if len(symsA) != len(symsB):
        raise AssertionError("Matrices use different numbers of symbols")

    canon = [sp.Symbol(f"x{i}") for i in range(len(symsA))]

    subsA = dict(zip(symsA, canon))
    subsB = dict(zip(symsB, canon))

    A2 = A.subs(subsA)
    B2 = B.subs(subsB)

    A2 = A2.applyfunc(sp.simplify)
    B2 = B2.applyfunc(sp.simplify)

    A2 = A2.applyfunc(lambda x: sp.nsimplify(x, rational=True))
    B2 = B2.applyfunc(lambda x: sp.nsimplify(x, rational=True))

    for a, b in zip(A2, B2):
        if sp.simplify(a - b) != 0:
            raise AssertionError(
                f"Matrices differ:\nA_simplified={A2}\nB_simplified={B2}\nOriginal A={A}\nOriginal B={B}"
            )


# ===============================================================
# 1. PLANE
# X(u,v) = (a u + b v, c u + d v, 0)
# ===============================================================


def test_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[1, 0], [0, 1]])
    assert_matrix_equal(M, expected)


def test_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[a**2 + c**2, a * b + c * d], [a * b + c * d, b**2 + d**2]])

    assert_matrix_equal(M, expected)


# ===============================================================
# 2. CYLINDER
# X(u,v) = (a cos u, a sin u, b v)
# ===============================================================


def test_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[R**2, 0], [0, 1]])
    assert_matrix_equal(M, expected)


def test_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[a**2, 0], [0, b**2]])
    assert_matrix_equal(M, expected)


# ===============================================================
# 3. SPHERE
# X(u,v) = (a sin u cos v, a sin u sin v, a cos u)
# ===============================================================


def test_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[R**2, 0], [0, R**2 * sp.sin(u) ** 2]])
    assert_matrix_equal(M, expected)


def test_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[a**2, 0], [0, a**2 * sp.sin(u) ** 2]])
    assert_matrix_equal(M, expected)


# ===============================================================
# 4. TORUS
# X(u,v) = ((a + b cos v) cos u, (a + b cos v) sin u, b sin v)
# ===============================================================


def test_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[(R + r * sp.cos(v)) ** 2, 0], [0, r**2]])
    assert_matrix_equal(M, expected)


def test_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[(a + b * sp.cos(v)) ** 2, 0], [0, b**2]])
    assert_matrix_equal(M, expected)


# ===============================================================
# 5. PARABOLOID
# X(u,v) = (u, v, a u² + b v²)
# ===============================================================


def test_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[1 + 4 * u**2, 4 * u * v], [4 * u * v, 1 + 4 * v**2]])
    assert_matrix_equal(M, expected)


def test_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [1 + (2 * a * u) ** 2, (2 * a * u) * (2 * b * v)],
            [(2 * a * u) * (2 * b * v), 1 + (2 * b * v) ** 2],
        ]
    )
    assert_matrix_equal(M, expected)


# ===============================================================
# 6. HYPERBOLIC PARABOLOID
# X(u,v) = (u, v, a u² − b v²)
# ===============================================================


def test_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[1 + 4 * u**2, -4 * u * v], [-4 * u * v, 1 + 4 * v**2]])
    assert_matrix_equal(M, expected)


def test_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [1 + (2 * a * u) ** 2, (2 * a * u) * (-2 * b * v)],
            [(2 * a * u) * (-2 * b * v), 1 + (2 * b * v) ** 2],
        ]
    )
    assert_matrix_equal(M, expected)


# ===============================================================
# 7. ELLIPTIC SURFACE (NON-ORTHOGONAL)
# X(u,v) = (u + v, a u − b v, a b u v)
# ===============================================================


def test_elliptic_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u + v, u - v, u * v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[2 + v**2, u * v], [u * v, 2 + u**2]])
    assert_matrix_equal(M, expected)


def test_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    X = parse_input(str(X), str([u, v]))
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    expected = sp.Matrix([[Xu.dot(Xu), Xu.dot(Xv)], [Xu.dot(Xv), Xv.dot(Xv)]])
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    assert_matrix_equal(M, expected)


# ===============================================================
# 8. CATENOID
# X(u,v) = (a cosh v cos u, a cosh v sin u, b v)
# ===============================================================


def test_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[sp.cosh(v) ** 2, 0], [0, sp.cosh(v) ** 2]])
    assert_matrix_equal(M, expected)


def test_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    E = a**2 * sp.cosh(v) ** 2
    F = 0
    G = a**2 * sp.sinh(v) ** 2 + b**2
    expected = sp.Matrix([[E, F], [F, G]])
    assert_matrix_equal(M, expected)


# ===============================================================
# 9. HELICOID
# X(u,v) = (a u cos v, a u sin v, b v)
# ===============================================================


def test_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[1, 0], [0, 1 + u**2]])
    assert_matrix_equal(M, expected)


def test_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    expected = sp.Matrix([[a**2, 0], [0, a**2 * u**2 + b**2]])
    assert_matrix_equal(M, expected)


# ===============================================================
# 10. NONLINEAR SURFACE
# X(u,v) = (a u² v, b u v², exp(a u + b v))
# ===============================================================


def test_weird_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    Xu = sp.Matrix([2 * u * v, v**2, sp.exp(u + v)])
    Xv = sp.Matrix([u**2, 2 * u * v, sp.exp(u + v)])
    expected = sp.Matrix([[Xu.dot(Xu), Xu.dot(Xv)], [Xu.dot(Xv), Xv.dot(Xv)]])
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    assert_matrix_equal(M, expected)


def test_weird_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    Xu = sp.Matrix([2 * a * u * v, b * v**2, a * sp.exp(a * u + b * v)])
    Xv = sp.Matrix([a * u**2, 2 * b * u * v, b * sp.exp(a * u + b * v)])
    expected = sp.Matrix([[Xu.dot(Xu), Xu.dot(Xv)], [Xu.dot(Xv), Xv.dot(Xv)]])
    M, _ = compute_first_fundamental_form(X, [u, v])
    M, _ = compute_first_fundamental_form(X, [u, v])
    assert_matrix_equal(M, expected)
