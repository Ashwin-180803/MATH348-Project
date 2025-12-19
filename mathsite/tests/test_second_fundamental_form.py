import sympy as sp
from calcapp.utils.parser import parse_input
import pytest
from calcapp.utils.functions import (
    compute_second_fundamental_form,
)


def assert_matrix_equal(A, B):
    assert A.shape == B.shape
    symsA = sorted(A.free_symbols, key=lambda s: s.name)
    symsB = sorted(B.free_symbols, key=lambda s: s.name)
    if len(symsA) != len(symsB):
        raise AssertionError("Matrices use different numbers of symbols")
    canon = [sp.Symbol(f"x{i}") for i in range(len(symsA))]
    subsA = dict(zip(symsA, canon))
    subsB = dict(zip(symsB, canon))
    A2 = A.subs(subsA).applyfunc(sp.simplify)
    B2 = B.subs(subsB).applyfunc(sp.simplify)
    A2 = A2.applyfunc(lambda x: sp.nsimplify(x, rational=True))
    B2 = B2.applyfunc(lambda x: sp.nsimplify(x, rational=True))
    for a, b in zip(A2, B2):
        if sp.simplify(a - b) != 0:
            raise AssertionError(
                f"Matrices differ:\nA_simplified={A2}\nB_simplified={B2}\nOriginal A={A}\nOriginal B={B}"
            )


# ===============================================================
# 1. PLANE
# ===============================================================
def test_plane_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[0, 0], [0, 0]])
    assert_matrix_equal(II, expected)


def test_plane_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[0, 0], [0, 0]])
    assert_matrix_equal(II, expected)


# ===============================================================
# 2. CYLINDER
# ===============================================================
def test_cylinder_numeric_second():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-R, 0], [0, 0]])
    assert_matrix_equal(II, expected)


def test_cylinder_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-a, 0], [0, 0]])
    assert_matrix_equal(II, expected)


# ===============================================================
# 3. SPHERE
# ===============================================================
def test_sphere_numeric_second():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-R, 0], [0, -R * sp.sin(u) ** 2]])
    II = II.applyfunc(sp.simplify)
    expected = expected.applyfunc(sp.simplify)
    assert_matrix_equal(II, expected)


def test_sphere_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-a, 0], [0, -a * sp.sin(u) ** 2]])
    assert_matrix_equal(II, expected)


# ===============================================================
# 4. TORUS
# X(u,v) = ((a + b cos v) cos u, (a + b cos v) sin u, b sin v)
# ===============================================================
def test_torus_numeric_second():
    u, v = sp.symbols("u v", real=True)
    R = 3.1
    r = 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-(R + r * sp.cos(v)) * sp.cos(v), 0], [0, -r]])
    II = II.applyfunc(sp.simplify)
    expected = expected.applyfunc(sp.simplify)
    assert_matrix_equal(II, expected)


def test_torus_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-(a + b * sp.cos(v)) * sp.cos(v), 0], [0, -b]])
    assert_matrix_equal(II, expected)


# ===============================================================
# 5. PARABOLOID
# X(u,v) = (u,v,a u^2 + b v^2)
# ===============================================================
def test_paraboloid_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [2 / sp.sqrt(4 * u**2 + 4 * v**2 + 1), 0],
            [0, 2 / sp.sqrt(4 * u**2 + 4 * v**2 + 1)],
        ]
    )
    assert_matrix_equal(II, expected)


def test_paraboloid_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [(2 * a) / sp.sqrt(4 * a**2 * u**2 + 4 * b**2 * v**2 + 1), 0],
            [0, (2 * b) / sp.sqrt(4 * a**2 * u**2 + 4 * b**2 * v**2 + 1)],
        ]
    )
    assert_matrix_equal(II, expected)


# ===============================================================
# 6. HYPERBOLIC PARABOLOID
# X(u,v) = (u,v,a u^2 - b v^2)
# ===============================================================
def test_hyperbolic_paraboloid_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [2 / sp.sqrt(4 * u**2 + 4 * v**2 + 1), 0],
            [0, -2 / sp.sqrt(4 * u**2 + 4 * v**2 + 1)],
        ]
    )
    assert_matrix_equal(II, expected)


def test_hyperbolic_paraboloid_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [(2 * a) / sp.sqrt(4 * a**2 * u**2 + 4 * b**2 * v**2 + 1), 0],
            [0, -(2 * b) / sp.sqrt(4 * a**2 * u**2 + 4 * b**2 * v**2 + 1)],
        ]
    )
    assert_matrix_equal(II, expected)


# ===============================================================
# 7. ELLIPTIC SURFACE
# X(u,v) = (u + v, a u − b v, a b u v)
# ===============================================================
def test_elliptic_surface_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u + v, u - v, u * v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])

    expected = sp.Matrix(
        [
            [0, -2 / sp.sqrt((u + v) ** 2 + (u - v) ** 2 + 4)],
            [-2 / sp.sqrt((u + v) ** 2 + (u - v) ** 2 + 4), 0],
        ]
    )
    assert_matrix_equal(II, expected)


def test_elliptic_surface_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u + v, a * u - b * v, a * b * u * v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])

    expected = sp.Matrix(
        [
            [
                0,
                -a
                * b
                * (a + b)
                / sp.sqrt(
                    (a * b * (a * u + b * v)) ** 2
                    + (a * b * (u - v)) ** 2
                    + (a + b) ** 2
                ),
            ],
            [
                -a
                * b
                * (a + b)
                / sp.sqrt(
                    (a * b * (a * u + b * v)) ** 2
                    + (a * b * (u - v)) ** 2
                    + (a + b) ** 2
                ),
                0,
            ],
        ]
    )
    assert_matrix_equal(II, expected)


# ===============================================================
# 8. CATENOID
# X(u,v) = (a cosh v cos u, a cosh v sin u, b v)
# ===============================================================
def test_catenoid_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[-1, 0], [0, 1]])

    assert_matrix_equal(II, expected)


def test_catenoid_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [-a * b * sp.cosh(v) / sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2), 0],
            [0, a * b * sp.cosh(v) / sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2)],
        ]
    )

    assert_matrix_equal(II, expected)


# ===============================================================
# 9. HELICOID
# X(u,v) = (a u cos v, a u sin v, b v)
# ===============================================================
def test_helicoid_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u * sp.cos(v), u * sp.sin(v), v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix([[0, -1 / sp.sqrt(u**2 + 1)], [-1 / sp.sqrt(u**2 + 1), 0]])

    assert_matrix_equal(II, expected)


def test_helicoid_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * u * sp.cos(v), a * u * sp.sin(v), b * v]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    expected = sp.Matrix(
        [
            [0, -(a * b) / sp.sqrt(a**2 * u**2 + b**2)],
            [-(a * b) / sp.sqrt(a**2 * u**2 + b**2), 0],
        ]
    )

    assert_matrix_equal(II, expected)


# ===============================================================
# 10. NONLINEAR SURFACE
# X(u,v) = (a u^2 v, b u v^2, exp(a u + b v))
# ===============================================================
def test_weird_surface_numeric_second():
    u, v = sp.symbols("u v", real=True)
    X = [u**2 * v, u * v**2, sp.exp(u + v)]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    rN_norm = sp.sqrt(
        ((v**2 - 2 * u * v) * sp.exp(u + v)) ** 2
        + ((u**2 - 2 * u * v) * sp.exp(u + v)) ** 2
        + (3 * u**2 * v**2) ** 2
    )

    expected = sp.Matrix(
        [
            [
                (
                    2 * v * (v**2 - 2 * u * v) * sp.exp(u + v)
                    + 0 * (u**2 - 2 * u * v) * sp.exp(u + v)
                    + sp.exp(u + v) * 3 * u**2 * v**2
                )
                / rN_norm,
                (
                    2 * u * (v**2 - 2 * u * v) * sp.exp(u + v)
                    + 2 * v * (u**2 - 2 * u * v) * sp.exp(u + v)
                    + sp.exp(u + v) * 3 * u**2 * v**2
                )
                / rN_norm,
            ],
            [
                (
                    2 * u * (v**2 - 2 * u * v) * sp.exp(u + v)
                    + 2 * v * (u**2 - 2 * u * v) * sp.exp(u + v)
                    + sp.exp(u + v) * 3 * u**2 * v**2
                )
                / rN_norm,
                (
                    0 * (v**2 - 2 * u * v) * sp.exp(u + v)
                    + 2 * u * (u**2 - 2 * u * v) * sp.exp(u + v)
                    + sp.exp(u + v) * 3 * u**2 * v**2
                )
                / rN_norm,
            ],
        ]
    )
    assert_matrix_equal(II, expected)


def test_weird_surface_symbolic_second():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]
    X = parse_input(str(X), str([u, v]))
    II, _ = compute_second_fundamental_form(X, [u, v])
    II, _ = compute_second_fundamental_form(X, [u, v])
    rN_norm = sp.sqrt(
        ((b**2 * v**2 - 2 * a * b * u * v) * sp.exp(a * u + b * v)) ** 2
        + ((a**2 * u**2 - 2 * a * b * u * v) * sp.exp(a * u + b * v)) ** 2
        + (3 * a * b * u**2 * v**2) ** 2
    )

    expected = sp.Matrix(
        [
            [
                (
                    2
                    * a
                    * v
                    * (b**2 * v**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + 0 * (a**2 * u**2 - 2 * a * b * u * v) * sp.exp(a * u + b * v)
                    + a**2 * sp.exp(a * u + b * v) * 3 * a * b * u**2 * v**2
                )
                / rN_norm,
                (
                    2
                    * a
                    * u
                    * (b**2 * v**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + 2
                    * b
                    * v
                    * (a**2 * u**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + a * b * sp.exp(a * u + b * v) * 3 * a * b * u**2 * v**2
                )
                / rN_norm,
            ],
            [
                (
                    2
                    * a
                    * u
                    * (b**2 * v**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + 2
                    * b
                    * v
                    * (a**2 * u**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + a * b * sp.exp(a * u + b * v) * 3 * a * b * u**2 * v**2
                )
                / rN_norm,
                (
                    0 * (b**2 * v**2 - 2 * a * b * u * v) * sp.exp(a * u + b * v)
                    + 2
                    * b
                    * u
                    * (a**2 * u**2 - 2 * a * b * u * v)
                    * sp.exp(a * u + b * v)
                    + b**2 * sp.exp(a * u + b * v) * 3 * a * b * u**2 * v**2
                )
                / rN_norm,
            ],
        ]
    )
    assert_matrix_equal(II, expected)
