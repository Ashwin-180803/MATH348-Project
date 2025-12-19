import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_gauss_equations


def assert_gauss_equalities(
    first_gauss_rhs,
    second_gauss_rhs,
    third_gauss_rhs,
    fourth_gauss_rhs,
    first_gauss_lhs,
    second_gauss_lhs,
    third_gauss_lhs,
    fourth_gauss_lhs,
):
    if sp.simplify(fourth_gauss_rhs - fourth_gauss_lhs) != 0:
        raise AssertionError(
            f"First gauss equality mismatch:\nLHS: {first_gauss_lhs}\nRHS: {first_gauss_rhs}"
        )
    if sp.simplify(second_gauss_rhs - second_gauss_lhs) != 0:
        raise AssertionError(
            f"Second gauss equality mismatch:\nLHS: {second_gauss_lhs}\nRHS: {second_gauss_rhs}"
        )
    if sp.simplify(third_gauss_rhs - third_gauss_lhs) != 0:
        raise AssertionError(
            f"Third gauss equality mismatch:\nLHS: {third_gauss_lhs}\nRHS: {third_gauss_rhs}"
        )
    if sp.simplify(fourth_gauss_rhs - fourth_gauss_lhs) != 0:
        raise AssertionError(
            f"Fourth gauss equality mismatch:\nLHS: {fourth_gauss_lhs}\nRHS: {fourth_gauss_rhs}"
        )


def assert_gauss_equal(
    computed_first_rhs,
    computed_second_rhs,
    computed_third_rhs,
    computed_fourth_rhs,
    expected_first_rhs,
    expected_second_rhs,
    expected_third_rhs,
    expected_fourth_rhs,
):
    computed_first_rhs = sp.sympify(computed_first_rhs)
    expected_first_rhs = sp.simplify(sp.sympify(expected_first_rhs))
    computed_second_rhs = sp.sympify(computed_second_rhs)
    expected_second_rhs = sp.simplify(sp.sympify(expected_second_rhs))
    computed_third_rhs = sp.sympify(computed_third_rhs)
    expected_third_rhs = sp.simplify(sp.sympify(expected_third_rhs))
    computed_fourth_rhs = sp.sympify(computed_fourth_rhs)
    expected_fourth_rhs = sp.simplify(sp.sympify(expected_fourth_rhs))

    first_eq_equal = sp.simplify(computed_first_rhs - expected_first_rhs) == 0
    second_eq_equal = sp.simplify(computed_second_rhs - expected_second_rhs) == 0
    third_eq_equal = sp.simplify(computed_third_rhs - expected_third_rhs) == 0
    fourth_eq_equal = sp.simplify(computed_fourth_rhs - expected_fourth_rhs) == 0

    if first_eq_equal and second_eq_equal and third_eq_equal and fourth_eq_equal:
        return
    else:
        raise AssertionError(
            f"Gauss equations do not match expected values:\n"
            f"First Gauss difference: {first_eq_equal}\n"
            f"Second Gauss difference: {second_eq_equal}\n"
            f"Third Gauss difference: {third_eq_equal}\n"
            f"Fourth Gauss difference: {fourth_eq_equal}"
        )


# ----------------------------
# 1. PLANE
# ----------------------------
def test_gauss_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 0
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = 0

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 0
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = 0

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


# ----------------------------
# 2. CYLINDER
# ----------------------------
def test_gauss_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2.7
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 0
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = 0

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 0
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = 0

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


# ----------------------------
# 3. SPHERE
# ----------------------------
def test_gauss_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = sp.nsimplify(1.8)
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 1
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = sp.sin(u) ** 2

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )

    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 1
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = sp.sin(u) ** 2

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )

    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


# ----------------------------
# 4. TORUS
# ----------------------------
def test_gauss_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R, r = sp.nsimplify(3.1), sp.nsimplify(0.8)
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = ((R + r * sp.cos(v)) * sp.cos(v)) / r
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = r * sp.cos(v) / (R + r * sp.cos(v))

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )

    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = ((a + b * sp.cos(v)) * sp.cos(v)) / b
    expected_second_rhs = 0
    expected_third_rhs = 0
    expected_fourth_rhs = b * sp.cos(v) / (a + b * sp.cos(v))

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


# ----------------------------
# 5. PARABOLOID
# ----------------------------
def test_gauss_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 + v**2]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = 4 * (1 + 4 * u**2) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_second_rhs = (16 * u * v) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_third_rhs = (16 * u * v) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_fourth_rhs = 4 * (1 + 4 * v**2) / ((1 + 4 * u**2 + 4 * v**2) ** 2)

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )

    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 + b * v**2]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])
    expected_first_rhs = (
        4 * a * b * (4 * a**2 * u**2 + 1) / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_second_rhs = (
        16 * a**2 * b**2 * u * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_third_rhs = (
        16 * a**2 * b**2 * u * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_fourth_rhs = (
        4 * a * b * (4 * b**2 * v**2 + 1) / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )

    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


# ----------------------------
# 6. HYPERBOLIC PARABOLOID
# ----------------------------
def test_gauss_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, u**2 - v**2]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = -4 * (1 + 4 * u**2) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_second_rhs = (16 * u * v) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_third_rhs = (16 * u * v) / ((1 + 4 * u**2 + 4 * v**2) ** 2)
    expected_fourth_rhs = -4 * (1 + 4 * v**2) / ((1 + 4 * u**2 + 4 * v**2) ** 2)

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )


def test_gauss_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [u, v, a * u**2 - b * v**2]
    X = parse_input(str(X), str([u, v]))
    (
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
        _,
        _,
    ) = compute_gauss_equations(X, [u, v])

    expected_first_rhs = (
        -4
        * a
        * b
        * (4 * a**2 * u**2 + 1)
        / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_second_rhs = (
        16 * a**2 * b**2 * u * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_third_rhs = (
        16 * a**2 * b**2 * u * v / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )
    expected_fourth_rhs = (
        -4
        * a
        * b
        * (4 * b**2 * v**2 + 1)
        / (4 * a**2 * u**2 + 4 * b**2 * v**2 + 1) ** 2
    )

    assert_gauss_equalities(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        first_lhs,
        second_lhs,
        third_lhs,
        fourth_lhs,
    )
    assert_gauss_equal(
        first_rhs,
        second_rhs,
        third_rhs,
        fourth_rhs,
        expected_first_rhs,
        expected_second_rhs,
        expected_third_rhs,
        expected_fourth_rhs,
    )
