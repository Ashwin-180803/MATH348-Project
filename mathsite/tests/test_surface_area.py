import sympy as sp
from calcapp.utils.parser import parse_input
import pytest
from calcapp.utils.functions import compute_first_fundamental_form, compute_surface_area


def assert_expr_equal(A, B):
    """
    Assert two sympy expressions are equal (symbolically).
    Uses simplify / nsimplify to handle rational floats, etc.
    """
    diff = sp.simplify(sp.nsimplify(A) - sp.nsimplify(B))
    if diff != 0:
        raise AssertionError(f"Expressions differ:\nA={A}\nB={B}\nA-B={diff}")

u1, u2, v1, v2 = sp.symbols("u1 u2 v1 v2", real=True)

# ===============================================================
# 1. PLANE
# X(u,v) = (a u + b v, c u + d v, 0)
# ===============================================================


def test_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = [u, v, 0]
    X = parse_input(str(X), str([u, v]))

    SA, _ = compute_surface_area(X, [u, v], [0, 1], [0, 1])
    SA, _ = compute_surface_area(X, [u, v], [0, 1], [0, 1])
    expected = sp.Integer(1)
    assert_expr_equal(SA, expected)


def test_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = [a * u + b * v, c * u + d * v, 0]
    X = parse_input(str(X), str([u, v]))

    M = compute_first_fundamental_form(X, [u, v])
    E = M[0, 0]
    F = M[0, 1]
    G = M[1, 1]
    integrand = sp.sqrt(E * G - F**2)

    expected = sp.integrate(integrand, (u, u1, u2), (v, v1, v2))
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])
    assert_expr_equal(SA, expected)


# ===============================================================
# 2. CYLINDER
# X(u,v) = (a cos u, a sin u, b v)
# ===============================================================


def test_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = sp.Rational(27, 10)  # 2.7
    h = sp.Integer(5)
    X = [R * sp.cos(u), R * sp.sin(u), v]
    X = parse_input(str(X), str([u, v]))

    SA, _ = compute_surface_area(X, [u, v], [0, 2 * sp.pi], [0, h])
    SA, _ = compute_surface_area(X, [u, v], [0, 2 * sp.pi], [0, h])
    expected = 2 * sp.pi * R * h
    assert_expr_equal(SA, expected)


def test_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(u), a * sp.sin(u), b * v]
    X = parse_input(str(X), str([u, v]))

    M = compute_first_fundamental_form(X, [u, v])
    integrand = sp.sqrt(M[0, 0] * M[1, 1] - M[0, 1] ** 2)
    expected = sp.integrate(integrand, (u, u1, u2), (v, v1, v2))
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])
    assert_expr_equal(SA, expected)


# ===============================================================
# 3. SPHERE
# X(u,v) = (a sin u cos v, a sin u sin v, a cos u)
# ===============================================================


def test_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = sp.Rational(18, 10)  # 1.8
    X = [R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))

    SA, _ = compute_surface_area(X, [u, v], [0, sp.pi], [0, 2 * sp.pi])
    SA, _ = compute_surface_area(X, [u, v], [0, sp.pi], [0, 2 * sp.pi])
    expected = 4 * sp.pi * R**2

    assert_expr_equal(SA, expected)


def test_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = [a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]
    X = parse_input(str(X), str([u, v]))
    SA, _ = compute_surface_area(X, [u, v], [0, sp.pi], [0, 2 * sp.pi])
    SA, _ = compute_surface_area(X, [u, v], [0, sp.pi], [0, 2 * sp.pi])

    expected = 4 * sp.pi * a**2

    assert_expr_equal(SA, expected)


# ===============================================================
# 4. TORUS
# X(u,v) = ((a + b cos v) cos u, (a + b cos v) sin u, b sin v)
# ===============================================================


def test_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = sp.Rational(31, 10)  # 3.1
    r = sp.Rational(8, 10)  # 0.8
    X = [
        (R + r * sp.cos(v)) * sp.cos(u),
        (R + r * sp.cos(v)) * sp.sin(u),
        r * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))

    SA, _ = compute_surface_area(X, [u, v], [0, 2 * sp.pi], [0, 2 * sp.pi])
    SA, _ = compute_surface_area(X, [u, v], [0, 2 * sp.pi], [0, 2 * sp.pi])
    expected = 4 * sp.pi**2 * R * r
    assert_expr_equal(SA, expected)


def test_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [
        (a + b * sp.cos(v)) * sp.cos(u),
        (a + b * sp.cos(v)) * sp.sin(u),
        b * sp.sin(v),
    ]
    X = parse_input(str(X), str([u, v]))

    expected = b * (u2 - u1) * (a * (v2 - v1) + b * (sp.sin(v2) - sp.sin(v1)))
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])
    SA, _ = compute_surface_area(X, [u, v], [u1, u2], [v1, v2])

    assert_expr_equal(SA, expected)
