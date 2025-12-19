import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_arc_length, true_simplify


def assert_arc_lengths_equal(v1, v2):
    """Compare two symbolic arc lengths robustly (scalar comparison)."""
    v1 = sp.sympify(v1)
    v2 = sp.sympify(v2)

    syms1 = sorted(v1.free_symbols, key=lambda s: s.name)
    syms2 = sorted(v2.free_symbols, key=lambda s: s.name)

    if len(syms1) != len(syms2):
        raise AssertionError("Different symbol counts in arc-length expressions")

    canon = [sp.Symbol(f"x{i}") for i in range(len(syms1))]

    v1 = v1.subs(dict(zip(syms1, canon)))
    v2 = v2.subs(dict(zip(syms2, canon)))

    v1 = sp.nsimplify(true_simplify(v1), rational=True)
    v2 = sp.nsimplify(true_simplify(v2), rational=True)

    diff = sp.simplify(v1 - v2)

    # For cases involving non-computable integrals
    if isinstance(v1, sp.Integral) and isinstance(v2, sp.Integral):
        # Extract integrand and limits
        f1, lims1 = v1.function, v1.limits
        f2, lims2 = v2.function, v2.limits

        if sp.simplify(f1 - f2) == 0:
            if lims1 == lims2:
                return
    elif diff != 0:
        raise AssertionError(f"Arc lengths differ: {v1} vs {v2}")


def expected_arc_length(X, t, bounds):
    X = [sp.sympify(expr) for expr in X]
    Xt = sp.Matrix([sp.diff(expr, t) for expr in X])
    integrand = sp.sqrt(sum(component**2 for component in Xt))
    t0, t1 = bounds
    return sp.simplify(sp.integrate(integrand, (t, t0, t1)))


# ===============================================================
# 1. LINE
# ===============================================================


def test_line_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([t, 2 * t, 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 3])
    ds, _ = compute_arc_length(X, t, [0, 3])
    expected = expected_arc_length([t, 2 * t, 0], t, [0, 3])

    assert_arc_lengths_equal(ds, expected)


def test_line_symbolic():
    t = sp.symbols("t", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = parse_input(str([a * t + b, c * t + d, 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * t + b, c * t + d, 0], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 2. CIRCLE
# ===============================================================


def test_circle_numeric():
    t = sp.symbols("t", real=True)
    R = 2.7
    X = parse_input(str([R * sp.cos(t), R * sp.sin(t), 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    expected = expected_arc_length([R * sp.cos(t), R * sp.sin(t), 0], t, [0, sp.pi])
    assert_arc_lengths_equal(ds, expected)


def test_circle_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = parse_input(str([a * sp.cos(t), a * sp.sin(t), 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * sp.cos(t), a * sp.sin(t), 0], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 3. SPHERE MERIDIAN
# ===============================================================


def test_sphere_meridian_numeric():
    t = sp.symbols("t", real=True)
    R = 1.8
    X = parse_input(str([R * sp.sin(t), 0, R * sp.cos(t)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, sp.pi / 2])
    ds, _ = compute_arc_length(X, t, [0, sp.pi / 2])
    expected = expected_arc_length([R * sp.sin(t), 0, R * sp.cos(t)], t, [0, sp.pi / 2])
    assert_arc_lengths_equal(ds, expected)


def test_sphere_meridian_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)

    X = parse_input(str([a * sp.sin(t), 0, a * sp.cos(t)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * sp.sin(t), 0, a * sp.cos(t)], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 4. TORUS CURVE
# ===============================================================


def test_torus_curve_numeric():
    t = sp.symbols("t", real=True)
    R = 3.1
    r = 0.8
    c = 1.2
    A = R + r * sp.cos(c)
    X = parse_input(str([A * sp.cos(t), A * sp.sin(t), r * sp.sin(c)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    expected = expected_arc_length(
        [A * sp.cos(t), A * sp.sin(t), r * sp.sin(c)], t, [0, sp.pi]
    )
    assert_arc_lengths_equal(ds, expected)


def test_torus_curve_symbolic():
    t = sp.symbols("t", real=True)
    a, b, c = sp.symbols("a b c", real=True, positive=True)
    A = a + b * sp.cos(c)
    X = parse_input(str([A * sp.cos(t), A * sp.sin(t), b * sp.sin(c)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length(
        [A * sp.cos(t), A * sp.sin(t), b * sp.sin(c)], t, [0, t]
    )
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 5. PARABOLA
# ===============================================================


def test_parabola_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([t, 0, t**2]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 1])
    ds, _ = compute_arc_length(X, t, [0, 1])
    expected = expected_arc_length([t, 0, t**2], t, [0, 1])
    assert_arc_lengths_equal(ds, expected)


def test_parabola_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = parse_input(str([t, 0, a * t**2]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([t, 0, a * t**2], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 6. HYPERBOLIC PARABOLA
# ===============================================================


def test_hyperbolic_parabola_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([t, 0, t**2 - t]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 2])
    ds, _ = compute_arc_length(X, t, [0, 2])
    expected = expected_arc_length([t, 0, t**2 - t], t, [0, 2])
    assert_arc_lengths_equal(ds, expected)


def test_hyperbolic_parabola_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([t, 0, a * t**2 - b * t]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([t, 0, a * t**2 - b * t], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 7. ELLIPTIC CURVE
# ===============================================================


def test_elliptic_curve_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([t, t, t**2]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 1])
    ds, _ = compute_arc_length(X, t, [0, 1])
    expected = expected_arc_length([t, t, t**2], t, [0, 1])
    assert_arc_lengths_equal(ds, expected)


def test_elliptic_curve_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([t, a * t, a * b * t**2]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([t, a * t, a * b * t**2], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 8. CATENARY
# ===============================================================


def test_catenary_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([sp.cosh(t), t, 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 1])
    ds, _ = compute_arc_length(X, t, [0, 1])
    expected = expected_arc_length([sp.cosh(t), t, 0], t, [0, 1])
    assert_arc_lengths_equal(ds, expected)


def test_catenary_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([a * sp.cosh(t), b * t, 0]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * sp.cosh(t), b * t], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 9. HELIX
# ===============================================================


def test_helix_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([sp.cos(t), sp.sin(t), 2 * t]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    ds, _ = compute_arc_length(X, t, [0, sp.pi])
    expected = expected_arc_length([sp.cos(t), sp.sin(t), 2 * t], t, [0, sp.pi])
    assert_arc_lengths_equal(ds, expected)


def test_helix_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([a * sp.cos(t), a * sp.sin(t), b * t]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * sp.cos(t), a * sp.sin(t), b * t], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 10. NONLINEAR CURVE
# ===============================================================


def test_weird_curve_numeric():
    t = sp.symbols("t", real=True)
    X = parse_input(str([t**2, t**3, sp.exp(t)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, 1])
    ds, _ = compute_arc_length(X, t, [0, 1])
    expected = expected_arc_length([t**2, t**3, sp.exp(t)], t, [0, 1])
    assert_arc_lengths_equal(ds, expected)


def test_weird_curve_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([a * t**2, b * t**3, sp.exp(a * t)]), str([t]))
    ds, _ = compute_arc_length(X, t, [0, t])
    ds, _ = compute_arc_length(X, t, [0, t])
    expected = expected_arc_length([a * t**2, b * t**3, sp.exp(a * t)], t, [0, t])
    assert_arc_lengths_equal(ds, expected)


# ===============================================================
# 11. DIFFERENT VARIABLES
# ===============================================================


def test_different_variables_numeric():
    m = sp.symbols("m", real=True)
    X = parse_input(str([m, m**2]), str([m]))
    ds, _ = compute_arc_length(X, m, [0, m])
    ds, _ = compute_arc_length(X, m, [0, m])
    expected = expected_arc_length([m, m**2], m, [0, m])

    assert_arc_lengths_equal(ds, expected)
