import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_shape_operator

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


# ======================================================================
# 1. PLANE
# X(u,v) = (u, v, 0)
# Shape operator = 0
# ======================================================================


def test_shape_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.zeros(2)
    assert_matrix_equal(S, expected)


def test_shape_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.zeros(2)
    assert_matrix_equal(S, expected)


# ======================================================================
# 2. CYLINDER
# X(u,v) = (a cos u, a sin u, b v)
# Principal curvatures: (1/a, 0)
# Shape operator = diag(1/a, 0)
# ======================================================================


def test_shape_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-1 / R, 0], [0, 0]])
    assert_matrix_equal(S, expected)


def test_shape_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-1 / a, 0], [0, 0]])
    assert_matrix_equal(S, expected)


# ======================================================================
# 3. SPHERE
# X(u,v) = (a sin u cos v, a sin u sin v, a cos u)
# All principal curvatures = 1/a
# Shape operator = (1/a) * I
# ======================================================================


def test_shape_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-1 / R, 0], [0, -1 / R]])
    assert_matrix_equal(S, expected)


def test_shape_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", positive=True, real=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-1 / a, 0], [0, -1 / a]])
    assert_matrix_equal(S, expected)


# ======================================================================
# 4. TORUS
# X(u,v) = ((a + b cos v)cos u, (a + b cos v) sin u, b sin v)
# Shape operator known:
#   S_11 = cos v / (a + b cos v)
#   S_22 = 1/b
# ======================================================================


def test_shape_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-sp.cos(v) / (R + r * sp.cos(v)), 0], [0, -1 / r]])
    assert_matrix_equal(S, expected)


def test_shape_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-sp.cos(v) / (a + b * sp.cos(v)), 0], [0, -1 / b]])
    assert_matrix_equal(S, expected)


# ======================================================================
# 5. PARABOLOID
# X(u,v) = (u, v, u² + v²)
# Known shape operator:
#   S = 1/(1 + 4u² + 4v²) * [[2, 0],[0, 2]]
# ======================================================================


def test_shape_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    factor = 2 / ((1 + 4 * u**2 + 4 * v**2) ** (3 / 2))
    expected = sp.Matrix(
        factor * sp.Matrix([[1 + 4 * v**2, -4 * u * v], [-4 * u * v, 1 + 4 * u**2]])
    )

    assert_matrix_equal(S, expected)


def test_shape_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    factor = 2 / ((1 + 4 * a**2 * u**2 + 4 * b**2 * v**2) ** (3 / 2))
    expected = sp.Matrix(
        factor
        * sp.Matrix(
            [
                [a * (1 + 4 * b**2 * v**2), -4 * a * b**2 * u * v],
                [-4 * a**2 * b * u * v, b * (1 + 4 * a**2 * u**2)],
            ]
        )
    )
    assert_matrix_equal(S, expected)


# ======================================================================
# 6. HYPERBOLIC PARABOLOID
# X(u,v) = (u, v, u² - v²)
# S = 1/(1 + 4u² + 4v²) * [[2, 0],[0, -2]]
# ======================================================================


def test_shape_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.simplify(
        sp.Matrix(
            [
                [1 + (2 * v) ** 2, 4 * u * v],
                [4 * u * v, 1 + (2 * u) ** 2],
            ]
        )
        @ sp.Matrix(
            [
                [2 * 1 / sp.sqrt(1 + 4 * u**2 + 4 * v**2), 0],
                [0, -2 * 1 / sp.sqrt(1 + 4 * u**2 + 4 * v**2)],
            ]
        )
        / ((1 + (2 * u) ** 2) * (1 + (2 * v) ** 2) - (-4 * u * v) ** 2)
    )
    assert_matrix_equal(S, expected)


def test_shape_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.simplify(
        sp.Matrix(
            [
                [1 + (2 * b * v) ** 2, 4 * a * b * u * v],
                [4 * a * b * u * v, 1 + (2 * a * u) ** 2],
            ]
        )
        @ sp.Matrix(
            [
                [2 * a / sp.sqrt(1 + 4 * a**2 * u**2 + 4 * b**2 * v**2), 0],
                [0, -2 * b / sp.sqrt(1 + 4 * a**2 * u**2 + 4 * b**2 * v**2)],
            ]
        )
        / ((1 + (2 * a * u) ** 2) * (1 + (2 * b * v) ** 2) - (-4 * a * b * u * v) ** 2)
    )
    assert_matrix_equal(S, expected)


# ======================================================================
# 7. ELLIPTIC TYPE SURFACE (NON-ORTHOGONAL)
# X(u,v) = (u+v, a u − b v, ab u v)
# ======================================================================


def test_shape_elliptic_surface_numeric():
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
    expected = I.inv() * II

    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    assert_matrix_equal(S, expected)


def test_shape_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    Xuu = sp.Matrix([0, 0, 0])
    Xuv = sp.Matrix([0, 0, a * b])
    Xvv = sp.Matrix([0, 0, 0])
    n = (Xu.cross(Xv)).normalized()
    L = Xuu.dot(n)
    M = Xuv.dot(n)
    N = Xvv.dot(n)
    E = Xu.dot(Xu)
    F = Xu.dot(Xv)
    G = Xv.dot(Xv)
    I = sp.Matrix([[E, F], [F, G]])
    II = sp.Matrix([[L, M], [M, N]])
    expected = I.inv() * II

    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    assert_matrix_equal(S, expected)


# ======================================================================
# 8. CATENOID
# X(u,v) = (a cosh v cos u, a cosh v sin u, b v)
# Principal curvatures: ±(1 / (a cosh² v))
# Shape operator = diag(1/(a cosh² v),  -1/(a cosh² v))
# ======================================================================


def test_shape_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[-1 / sp.cosh(v) ** 2, 0], [0, 1 / sp.cosh(v) ** 2]])
    assert_matrix_equal(S, expected)


def test_shape_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix(
        [
            [-b / (a * sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2) * sp.cosh(v)), 0],
            [0, a * b * sp.cosh(v) / (a**2 * sp.sinh(v) ** 2 + b**2) ** (3 / 2)],
        ]
    )
    assert_matrix_equal(S, expected)


# ======================================================================
# 9. HELICOID
# X(u,v) = (a u cos v, a u sin v, b v)
# Principal curvatures: ±(b / (a² u² + b²))
# ======================================================================


def test_shape_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix([[0, -1 / sp.sqrt(1 + u**2)], [-1 / (1 + u**2) ** (3 / 2), 0]])
    assert_matrix_equal(S, expected)


def test_shape_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    expected = sp.Matrix(
        [
            [0, -b / (a * (sp.sqrt(b**2 + a**2 * u**2)))],
            [-(a * b) / (b**2 + a**2 * u**2) ** (3 / 2), 0],
        ]
    )

    assert_matrix_equal(S, expected)


# ======================================================================
# 10. NONLINEAR SURFACE
# X(u,v) = (a u² v, b u v², exp(a u + b v))
# ======================================================================


def test_shape_weird_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    Xu = sp.Matrix([2 * u * v, v**2, sp.exp(u + v)])
    Xv = sp.Matrix([u**2, 2 * u * v, sp.exp(u + v)])
    Xuu = sp.Matrix([2 * v, 0, sp.exp(u + v)])
    Xuv = sp.Matrix([2 * u, 2 * v, sp.exp(u + v)])
    Xvv = sp.Matrix([0, 2 * u, sp.exp(u + v)])
    n = (Xu.cross(Xv)).normalized()
    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)
    expected = sp.Matrix([[E, F], [F, G]]).inv() * sp.Matrix([[L, M], [M, N]])

    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    assert_matrix_equal(S, expected)


def test_shape_weird_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", positive=True, real=True)
    Xu = sp.Matrix([2 * a * u * v, b * v**2, a * sp.exp(a * u + b * v)])
    Xv = sp.Matrix([a * u**2, 2 * b * u * v, b * sp.exp(a * u + b * v)])
    Xuu = sp.Matrix([2 * a * v, 0, a**2 * sp.exp(a * u + b * v)])
    Xuv = sp.Matrix([2 * a * u, 2 * b * v, a * b * sp.exp(a * u + b * v)])
    Xvv = sp.Matrix([0, 2 * b * u, b**2 * sp.exp(a * u + b * v)])
    n = (Xu.cross(Xv)).normalized()

    E, F, G = Xu.dot(Xu), Xu.dot(Xv), Xv.dot(Xv)
    L, M, N = Xuu.dot(n), Xuv.dot(n), Xvv.dot(n)

    expected = sp.Matrix([[E, F], [F, G]]).inv() * sp.Matrix([[L, M], [M, N]])

    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    S, _ = compute_shape_operator(X, [u, v])
    S, _ = compute_shape_operator(X, [u, v])
    assert_matrix_equal(S, expected)
