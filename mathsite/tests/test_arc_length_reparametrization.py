import sympy as sp
import pytest

from calcapp.utils.parser import parse_input
from calcapp.utils.functions import (
    compute_arc_length_reparametrization,
    true_simplify,
)


def assert_arc_lengths_equal(result, expected):
    """
    Assert that two parametric curves result(s) and expected(s)
    have identical arc-length speed functions with respect to s.
    """
    s = sp.symbols("s", real=True)

    def speed(curve):
        deriv = [sp.diff(comp, s) for comp in curve]
        return sp.simplify(sp.sqrt(sum(d**2 for d in deriv)))

    L_res = speed(result)
    L_exp = speed(expected)

    if not sp.simplify(L_res - L_exp) == 0:
        raise AssertionError(
            f"Arc-length mismatch:\n"
            f"  result speed   = {L_res}\n"
            f"  expected speed = {L_exp}"
        )
    return True


# ======================================================================
# 1. STRAIGHT LINE
# X(t) = (a t, b t, c t)
# ======================================================================


def test_line_numeric():
    t = sp.symbols("t", real=True)
    X = [3 * t, 4 * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)

    s = result[0].free_symbols.pop()
    expected = [3 * s / 5, 4 * s / 5, 0]

    assert_arc_lengths_equal(result, expected)


def test_line_symbolic():
    t = sp.symbols("t", real=True)
    a, b, c = sp.symbols("a b c", real=True, positive=True)
    X = [a * t, b * t, c * t]
    bounds = [0, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [
        a * s / sp.sqrt(a**2 + b**2 + c**2),
        b * s / sp.sqrt(a**2 + b**2 + c**2),
        c * s / sp.sqrt(a**2 + b**2 + c**2),
    ]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 2. CIRCLE ARC
# X(t) = (R cos t, R sin t, 0)
# ======================================================================


def test_circle_numeric():
    t = sp.symbols("t", real=True)
    R = 2.3
    X = [R * sp.cos(t), R * sp.sin(t), 0]
    bounds = [1, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [R * sp.cos(s / R + 1), R * sp.sin(s / R + 1), 0]

    assert_arc_lengths_equal(result, expected)


def test_circle_symbolic():
    t = sp.symbols("t", real=True)
    R = sp.symbols("R", real=True, positive=True)
    X = [R * sp.cos(t), R * sp.sin(t), 0]
    bounds = [1, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [R * sp.cos(s / R + 1), R * sp.sin(s / R + 1), 0]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 3. HELIX
# X(t) = (a cos t, a sin t, b t)
# ======================================================================


def test_helix_numeric():
    t = sp.symbols("t", real=True)
    a = 1.1
    b = 0.9
    X = [a * sp.cos(t), a * sp.sin(t), b * t]
    bounds = [-2, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [
        11 * sp.cos(5 * sp.sqrt(202) * s / 101 - 2) / 10,
        11 * sp.sin(5 * sp.sqrt(202) * s / 101 - 2) / 10,
        9 * sp.sqrt(202) * s / 202 - 9 / 5,
    ]

    assert_arc_lengths_equal(result, expected)


def test_helix_symbolic():
    t = sp.symbols("t", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = [a * sp.cos(t), a * sp.sin(t), b * t]
    bounds = [-2, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    L = sp.sqrt(a**2 + b**2)
    expected = [a * sp.cos((s / L) - 2), a * sp.sin((s / L) - 2), b * ((s / L) - 2)]
    expected = [a * sp.cos((s / L) - 2), a * sp.sin((s / L) - 2), b * ((s / L) - 2)]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 4. PARABOLA
# X(t) = (t, t^2, 0)
# ======================================================================


def test_parabola_numeric():
    t = sp.symbols("t", real=True)
    X = [t, sp.pi * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [s / sp.sqrt(1 + sp.pi**2), sp.pi * s / sp.sqrt(1 + sp.pi**2), 0]

    assert_arc_lengths_equal(result, expected)


def test_parabola_symbolic():
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = [t, a * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [s / sp.sqrt(1 + a**2), a * s / sp.sqrt(1 + a**2), 0]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 5. EXPONENTIAL CURVE
# X(t) = (exp(a t), exp(b t), 0)
# ======================================================================


def test_exponential_numeric():
    t = sp.symbols("t", real=True)
    X = [sp.exp(2 * t), sp.exp(3 * t), 0]
    bounds = [sp.log(2), t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [
        (729 * s**2 + 4320 * sp.sqrt(1, 1) * sp.sqrt(10) * s + 64000) ** (1 / 3) / 9
        - 4 / 9,
        ((729 * s**2 + 4320 * sp.sqrt(10) * s + 64000) ** (1 / 3) - 4) ** (3 / 2) / 27,
        0,
    ]

    assert_arc_lengths_equal(result, expected)


def test_exponential_symbolic():
    t = sp.symbols("t", real=True)
    a, b, c = sp.symbols("a b c", real=True, positive=True)
    X = [a * sp.exp(c * t), b * sp.exp(c * t), 0]
    bounds = [sp.log(2), t]

    Xp = parse_input(str(X), str([t]))
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    result, _ = compute_arc_length_reparametrization(Xp, t, bounds)
    s = result[0].free_symbols.pop()
    expected = [
        a * (s / sp.sqrt(a**2 + b**2) + 2**c),
        b * (s / sp.sqrt(a**2 + b**2) + 2**c),
        0,
    ]

    assert_arc_lengths_equal(result, expected)
