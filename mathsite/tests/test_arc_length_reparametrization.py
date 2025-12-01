import sympy as sp
import pytest

from calcapp.utils.parser import parse_input
from calcapp.utils.functions import (
    compute_arc_length_reparametrization,
    true_simplify,
)


# --------------------------------------------------------------
# Helper: arc length scalar comparison
# --------------------------------------------------------------
def assert_arc_lengths_equal(v1, v2):
    """
    Compare two lists of symbolic expressions robustly (elementwise).
    Canonicalizes symbols to compare expressions with different names.
    """
    if len(v1) != len(v2):
        raise AssertionError("Different lengths")

    # Flatten expressions
    v1_exprs = [sp.sympify(expr) for expr in v1]
    v2_exprs = [sp.sympify(expr) for expr in v2]

    # Collect all free symbols from both lists
    syms1 = sorted(
        {s for expr in v1_exprs for s in expr.free_symbols}, key=lambda s: s.name
    )
    syms2 = sorted(
        {s for expr in v2_exprs for s in expr.free_symbols}, key=lambda s: s.name
    )

    if len(syms1) != len(syms2):
        raise AssertionError("Different number of symbols")

    # Create canonical symbols
    canon = [sp.Symbol(f"x{i}") for i in range(len(syms1))]
    subs1 = dict(zip(syms1, canon))
    subs2 = dict(zip(syms2, canon))

    # Substitute and simplify elementwise
    for a, b in zip(v1_exprs, v2_exprs):
        a_sub = sp.nsimplify(true_simplify(a.subs(subs1)), rational=False)
        b_sub = sp.nsimplify(true_simplify(b.subs(subs2)), rational=False)

        diff = sp.simplify(a_sub - b_sub)

        # Handle identical Integrals
        if isinstance(a_sub, sp.Integral) and isinstance(b_sub, sp.Integral):
            if a_sub.function == b_sub.function and a_sub.limits == b_sub.limits:
                continue
            else:
                raise AssertionError(f"Arc lengths differ: {a_sub} vs {b_sub}")
        elif diff != 0:
            raise AssertionError(f"Arc lengths differ: {a_sub} vs {b_sub}")


# ======================================================================
# 1. STRAIGHT LINE
# X(t) = (a t, b t, c t)
# ======================================================================


def test_line_numeric():
    t, s = sp.symbols("t s", real=True)
    X = [3 * t, 4 * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [3 * s / 5, 4 * s / 5, 0]
    assert_arc_lengths_equal(result, expected)


def test_line_symbolic():
    t, a, b, c, s = sp.symbols("t a b c s", real=True)
    X = [a * t, b * t, c * t]
    bounds = [0, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [
        a * s / sp.sqrt(a**2 + b**2 + c**2),
        b * s / sp.sqrt(a**2 + b**2 + c**2),
        c * s / sp.sqrt(a**2 + b**2 + c**2),
    ]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 2. CIRCLE ARC
# X(t) = (R cos t, R sin t, 0)
# bounds shifted (1 to t)
# ======================================================================


def test_circle_numeric():
    t, s = sp.symbols("t s", real=True)
    R = 2.3
    X = [R * sp.cos(t), R * sp.sin(t), 0]
    bounds = [1, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [R * sp.cos(s / R + 1), R * sp.sin(s / R + 1), 0]

    assert_arc_lengths_equal(result, expected)


def test_circle_symbolic():
    t, R, s = sp.symbols("t R s", real=True, positive=True)
    X = [R * sp.cos(t), R * sp.sin(t), 0]
    bounds = [1, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [R * sp.cos(s / R + 1), R * sp.sin(s / R + 1), 0]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 3. HELIX
# X(t) = (a cos t, a sin t, b t)
# bounds [-2, t]
# ======================================================================


def test_helix_numeric():
    t, s = sp.symbols("t s", real=True)
    a = 1.1
    b = 0.9
    X = [a * sp.cos(t), a * sp.sin(t), b * t]
    bounds = [-2, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)

    expected = [
        11 * sp.cos(5 * sp.sqrt(202) * s / 101 - 2) / 10,
        11 * sp.sin(5 * sp.sqrt(202) * s / 101 - 2) / 10,
        9 * sp.sqrt(202) * s / 202 - 9 / 5,
    ]

    assert_arc_lengths_equal(result, expected)


def test_helix_symbolic():
    t, a, b, s = sp.symbols("t a b s", real=True)
    X = [a * sp.cos(t), a * sp.sin(t), b * t]
    bounds = [-2, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    L = sp.sqrt(a**2 + b**2)
    expected = [a * sp.cos(s / L - 2), a * sp.sin(s / L - 2), b * (s / L - 2)]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 4. PARABOLA
# X(t) = (t, t^2, 0)
# bounds  [0, 3]
# ======================================================================


def test_parabola_numeric():
    t, s = sp.symbols("t s", real=True)
    X = [t, sp.pi * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [s / sp.sqrt(1 + sp.pi**2), sp.pi * s / sp.sqrt(1 + sp.pi**2), 0]

    assert_arc_lengths_equal(result, expected)


def test_parabola_symbolic():
    t, a, s = sp.symbols("t a s", real=True)
    X = [t, a * t, 0]
    bounds = [0, t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [s / sp.sqrt(1 + a**2), a * s / sp.sqrt(1 + a**2), 0]

    assert_arc_lengths_equal(result, expected)


# ======================================================================
# 5. EXPONENTIAL CURVE
# X(t) = (exp(a t), exp(b t), 0)
# bounds [2, 5]
# ======================================================================


def test_exponential_numeric():
    t, s = sp.symbols("t s", real=True)
    X = [sp.exp(2 * t), sp.exp(3 * t), 0]
    bounds = [sp.log(2), t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [
        (729 * s**2 + 4320 * sp.sqrt(1, 1) * sp.sqrt(10) * s + 64000) ** (1 / 3) / 9
        - 4 / 9,
        ((729 * s**2 + 4320 * sp.sqrt(10) * s + 64000) ** (1 / 3) - 4) ** (3 / 2) / 27,
        0,
    ]

    assert_arc_lengths_equal(result, expected)


def test_exponential_symbolic():
    t, a, b, c, s = sp.symbols("t a b c s", real=True)
    X = [a * sp.exp(c * t), b * sp.exp(c * t), 0]
    bounds = [sp.log(2), t]

    Xp = parse_input(str(X))
    result = compute_arc_length_reparametrization(Xp, t, bounds)
    expected = [
        a * (s / sp.sqrt(a**2 + b**2) + 2**c),
        b * (s / sp.sqrt(a**2 + b**2) + 2**c),
        0,
    ]

    assert_arc_lengths_equal(result, expected)
