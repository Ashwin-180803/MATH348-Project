import sympy as sp
import pytest
from calcapp.utils.parser import parse_input
from calcapp.utils.functions import compute_unit_normal_vector


def assert_unit_normals_equal(n1, n2):
    n1 = sp.Matrix(n1)
    n2 = sp.Matrix(n2)

    if n1.shape != n2.shape:
        raise AssertionError(f"Shape mismatch: {n1.shape} vs {n2.shape}")

    n1 = n1 / sp.sqrt(sum(c**2 for c in n1))
    n2 = n2 / sp.sqrt(sum(c**2 for c in n2))

    syms1 = sorted(n1.free_symbols, key=lambda s: s.name)
    syms2 = sorted(n2.free_symbols, key=lambda s: s.name)
    if len(syms1) != len(syms2):
        raise AssertionError("Different symbol counts")

    canon = [sp.Symbol(f"x{i}") for i in range(len(syms1))]
    n1 = n1.subs(dict(zip(syms1, canon)))
    n2 = n2.subs(dict(zip(syms2, canon)))

    n1 = n1.applyfunc(sp.simplify).applyfunc(lambda x: sp.nsimplify(x, rational=True))
    n2 = n2.applyfunc(sp.simplify).applyfunc(lambda x: sp.nsimplify(x, rational=True))

    same = all(sp.simplify(a - b) == 0 for a, b in zip(n1, n2))
    opp = all(sp.simplify(a + b) == 0 for a, b in zip(n1, n2))

    if not (same or opp):
        raise AssertionError(f"Normals differ (even up to sign): n1={n1}, n2={n2}")


# ====================================================================
# 1. PLANE
# X(u,v) = (u, v, 0)
# ====================================================================


def test_plane_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input("[u, v, 0]", str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([0, 0, 1])
    assert_unit_normals_equal(N, expected)


def test_plane_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b, c, d = sp.symbols("a b c d", real=True, positive=True)
    X = parse_input(str([a * u + b * v, c * u + d * v, 0]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])

    expected = sp.Matrix([0, 0, 1])

    assert_unit_normals_equal(N, expected)


# ====================================================================
# 2. CYLINDER
# X(u,v) = (R cos u, R sin u, v)
# ====================================================================


def test_cylinder_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3
    X = parse_input(str([R * sp.cos(u), R * sp.sin(u), v]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.cos(u), sp.sin(u), 0])
    assert_unit_normals_equal(N, expected)


def test_cylinder_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([a * sp.cos(u), a * sp.sin(u), b * v]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.cos(u), sp.sin(u), 0])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 3. SPHERE
# X(u,v) = (a sin u cos v, a sin u sin v, a cos u)
# ====================================================================


def test_sphere_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 2
    X = parse_input(
        str([R * sp.sin(u) * sp.cos(v), R * sp.sin(u) * sp.sin(v), R * sp.cos(u)]),
        str([u, v]),
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.sin(u) * sp.cos(v), sp.sin(u) * sp.sin(v), sp.cos(u)])
    assert_unit_normals_equal(N, expected)


def test_sphere_symbolic():
    u, v = sp.symbols("u v", real=True)
    a = sp.symbols("a", real=True, positive=True)
    X = parse_input(
        str([a * sp.sin(u) * sp.cos(v), a * sp.sin(u) * sp.sin(v), a * sp.cos(u)]),
        str([u, v]),
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.sin(u) * sp.cos(v), sp.sin(u) * sp.sin(v), sp.cos(u)])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 4. TORUS
# X(u,v) = ((a + b cos v) cos u, ...)
# ====================================================================


def test_torus_numeric():
    u, v = sp.symbols("u v", real=True)
    R = 3
    r = 1
    X = parse_input(
        str(
            [
                (R + r * sp.cos(v)) * sp.cos(u),
                (R + r * sp.cos(v)) * sp.sin(u),
                r * sp.sin(v),
            ]
        ),
        str([u, v]),
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.cos(u) * sp.cos(v), sp.sin(u) * sp.cos(v), sp.sin(v)])
    assert_unit_normals_equal(N, expected)


def test_torus_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(
        str(
            [
                (a + b * sp.cos(v)) * sp.cos(u),
                (a + b * sp.cos(v)) * sp.sin(u),
                b * sp.sin(v),
            ]
        ),
        str([u, v]),
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([sp.cos(u) * sp.cos(v), sp.sin(u) * sp.cos(v), sp.sin(v)])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 5. PARABOLOID
# X(u,v) = (u, v, u² + v²)
# ====================================================================


def test_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(str([u, v, u**2 + v**2]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([-2 * u, -2 * v, 1])
    assert_unit_normals_equal(N, expected)


def test_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([u, v, a * u**2 + b * v**2]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([-2 * a * u, -2 * b * v, 1])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 6. HYPERBOLIC PARABOLOID
# X(u,v) = (u, v, u² − v²)
# ====================================================================


def test_hyperbolic_paraboloid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(str([u, v, u**2 - v**2]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([-2 * u, 2 * v, 1])
    assert_unit_normals_equal(N, expected)


def test_hyperbolic_paraboloid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([u, v, a * u**2 - b * v**2]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    expected = sp.Matrix([-2 * a * u, 2 * b * v, 1])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 7. ELLIPTIC SURFACE
# X(u,v) = (u+v, a u − b v, a b u v)
# ====================================================================


def test_elliptic_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(str([u + v, u - v, u * v]), str([u, v]))
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    Xu = sp.Matrix([1, 1, v])
    Xv = sp.Matrix([1, -1, u])
    expected = Xu.cross(Xv)
    assert_unit_normals_equal(N, expected)


def test_elliptic_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([u + v, a * u - b * v, a * b * u * v]), str([u, v]))
    Xu = sp.Matrix([1, a, a * b * v])
    Xv = sp.Matrix([1, -b, a * b * u])
    expected = Xu.cross(Xv)
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    assert_unit_normals_equal(N, expected)


# ====================================================================
# 8. CATENOID
# X(u,v) = (a cosh v cos u, ..., b v)
# ====================================================================


def test_catenoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(
        str([sp.cosh(v) * sp.cos(u), sp.cosh(v) * sp.sin(u), v]), str([u, v])
    )
    expected = sp.Matrix(
        [[sp.cos(u) / sp.cosh(v)], [sp.sin(u) / sp.cosh(v)], [-sp.tanh(v)]]
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    assert_unit_normals_equal(N, expected)


def test_catenoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(
        str([a * sp.cosh(v) * sp.cos(u), a * sp.cosh(v) * sp.sin(u), b * v]),
        str([u, v]),
    )
    expected = sp.Matrix(
        [
            [b * sp.cos(u) / sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2)],
            [b * sp.sin(u) / sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2)],
            [-a * sp.sinh(v) / sp.sqrt(a**2 * sp.sinh(v) ** 2 + b**2)],
        ]
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])

    assert_unit_normals_equal(N, expected)


# ====================================================================
# 9. HELICOID
# X(u,v) = (a u cos v, a u sin v, b v)
# ====================================================================


def test_helicoid_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(str([u * sp.cos(v), u * sp.sin(v), v]), str([u, v]))
    Xu = sp.Matrix([sp.cos(v), sp.sin(v), 0])
    Xv = sp.Matrix([-u * sp.sin(v), u * sp.cos(v), 1])
    expected = Xu.cross(Xv)
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    assert_unit_normals_equal(N, expected)


def test_helicoid_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(str([a * u * sp.cos(v), a * u * sp.sin(v), b * v]), str([u, v]))
    expected = sp.Matrix(
        [
            [b * sp.sin(v) / sp.sqrt(a**2 * u**2 + b**2)],
            [-b * sp.cos(v) / sp.sqrt(a**2 * u**2 + b**2)],
            [a * u / sp.sqrt(a**2 * u**2 + b**2)],
        ]
    )
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])

    assert_unit_normals_equal(N, expected)


# ====================================================================
# 10. WEIRD SURFACE
# X(u,v) = (a u² v, b u v², exp(a u + b v))
# ====================================================================


def test_weird_surface_numeric():
    u, v = sp.symbols("u v", real=True)
    X = parse_input(str([u**2 * v, u * v**2, sp.exp(u + v)]), str([u, v]))
    Xu = sp.Matrix([2 * u * v, v**2, sp.exp(u + v)])
    Xv = sp.Matrix([u**2, 2 * u * v, sp.exp(u + v)])
    expected = Xu.cross(Xv)
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    assert_unit_normals_equal(N, expected)


def test_weird_surface_symbolic():
    u, v = sp.symbols("u v", real=True)
    a, b = sp.symbols("a b", real=True, positive=True)
    X = parse_input(
        str([a * u**2 * v, b * u * v**2, sp.exp(a * u + b * v)]), str([u, v])
    )
    Xu = sp.Matrix([2 * a * u * v, b * v**2, a * sp.exp(a * u + b * v)])
    Xv = sp.Matrix([a * u**2, 2 * b * u * v, b * sp.exp(a * u + b * v)])
    expected = Xu.cross(Xv)
    N, _ = compute_unit_normal_vector(X, [u, v])
    N, _ = compute_unit_normal_vector(X, [u, v])
    assert_unit_normals_equal(N, expected)
