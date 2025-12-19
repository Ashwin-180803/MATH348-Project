import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_codazzi_equations


def assert_codazzi_equalities(
    first_codazzi_rhs, second_codazzi_rhs, first_codazzi_lhs, second_codazzi_lhs
):
    if not sp.Eq(first_codazzi_lhs, first_codazzi_rhs):
        raise AssertionError(
            f"First Codazzi equality mismatch:\nLHS: {first_codazzi_lhs}\nRHS: {first_codazzi_rhs}"
        )
    if not sp.Eq(second_codazzi_lhs, second_codazzi_rhs):
        raise AssertionError(
            f"Second Codazzi equality mismatch:\nLHS: {second_codazzi_lhs}\nRHS: {second_codazzi_rhs}"
        )


def assert_codazzi_equal(
    computed_first_rhs, computed_second_rhs, expected_first_rhs, expected_second_rhs
):
    computed_first_rhs = sp.sympify(computed_first_rhs)
    expected_first_rhs = sp.simplify(sp.sympify(expected_first_rhs))
    computed_second_rhs = sp.sympify(computed_second_rhs)
    expected_second_rhs = sp.simplify(sp.sympify(expected_second_rhs))

    first_eq_equal = computed_first_rhs.equals(expected_first_rhs)
    second_eq_equal = computed_second_rhs.equals(expected_second_rhs)

    if first_eq_equal and second_eq_equal:
        return
    else:
        raise AssertionError(
            f"Codazzi equations do not match expected values:\n"
            f"First Codazzi difference: {first_eq_equal}\n"
            f"Second Codazzi difference: {second_eq_equal}"
        )


# ----------------------------
# 1. PLANE
# ----------------------------
def test_codazzi_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


# ----------------------------
# 2. CYLINDER
# ----------------------------
def test_codazzi_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))

    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


# ----------------------------
# 3. SPHERE
# ----------------------------
def test_codazzi_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = R * 2 * sp.sin(u) * sp.cos(u)

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = 0
    expected_second_rhs = a * 2 * sp.sin(u) * sp.cos(u)

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


# ----------------------------
# 4. TORUS
# ----------------------------
def test_codazzi_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R, r = 3.1, 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = (R + 2 * r * sp.cos(v)) * sp.sin(v)
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = (a + 2 * b * sp.cos(v)) * sp.sin(v)
    expected_second_rhs = 0

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


# ----------------------------
# 5. PARABOLOID
# ----------------------------
def test_codazzi_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = -8 * v / (4 * u**2 + 4 * v**2 + 1) ** (3 / 2)
    expected_second_rhs = 8 * u / (4 * u**2 + 4 * v**2 + 1) ** (3 / 2)

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = (
        -8 * a * b**2 * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** (3 / 2)
    )
    expected_second_rhs = (
        8 * a**2 * b * u / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** (3 / 2)
    )

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


# ----------------------------
# 6. HYPERBOLIC PARABOLOID
# ----------------------------
def test_codazzi_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = -8 * v / (4 * u**2 + 4 * v**2 + 1) ** (3 / 2)
    expected_second_rhs = -8 * u / (4 * u**2 + 4 * v**2 + 1) ** (3 / 2)

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)


def test_codazzi_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )
    first_rhs, second_rhs, first_lhs, second_lhs, _ = compute_codazzi_equations(
        X, [u, v]
    )

    expected_first_rhs = (
        -8 * a * b**2 * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** (3 / 2)
    )
    expected_second_rhs = (
        -8 * a**2 * b * u / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** (3 / 2)
    )

    assert_codazzi_equalities(first_rhs, second_rhs, first_lhs, second_lhs)
    assert_codazzi_equal(first_rhs, second_rhs, expected_first_rhs, expected_second_rhs)
