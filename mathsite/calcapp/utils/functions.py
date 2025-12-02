import numpy as np
import sympy as sp
from sympy import Matrix, Q, integrate, sympify, Expr, assuming
from sympy.physics.units.quantities import Quantity
from func_timeout import func_timeout, FunctionTimedOut
import math
import os


# Curves


def numeric_curve_positions(curve_name, params, t0, t1, n):
    t = np.linspace(t0, t1, n)

    if curve_name == "line":
        P = np.array(
            [
                params.get("x0", 0.0),
                params.get("y0", 0.0),
                params.get("z0", 0.0),
            ],
            dtype=float,
        )
        Q = np.array(
            [
                params.get("x1", 1.0),
                params.get("y1", 0.0),
                params.get("z1", 0.0),
            ],
            dtype=float,
        )
        R = P[np.newaxis, :] + np.outer(t, (Q - P))

    elif curve_name == "circle":
        a = params.get("a", 1.0)
        R = np.column_stack((a * np.cos(t), a * np.sin(t), np.zeros_like(t)))

    elif curve_name == "ellipse":
        a = params.get("a", 2.0)
        b = params.get("b", 1.0)
        R = np.column_stack((a * np.cos(t), b * np.sin(t), np.zeros_like(t)))

    elif curve_name == "helix":
        a = params.get("a", 1.0)
        b = params.get("b", 0.2)
        R = np.column_stack((a * np.cos(t), a * np.sin(t), b * t))

    elif curve_name == "cycloid":
        a = params.get("a", 1.0)
        R = np.column_stack(
            (a * (t - np.sin(t)), a * (1 - np.cos(t)), np.zeros_like(t))
        )

    elif curve_name == "twisted_cubic":
        R = np.column_stack((t, t**2, t**3))

    elif curve_name == "catenary":
        C = params.get("C", 1.0)
        R = np.column_stack((t, C * np.cosh(t / C), np.zeros_like(t)))

    elif curve_name == "hyperbola":
        R = np.column_stack((np.cosh(t), np.sinh(t), np.zeros_like(t)))

    elif curve_name == "tractrix":
        R = np.column_stack((t - np.tanh(t), 1.0 / np.cosh(t), np.zeros_like(t)))

    elif curve_name == "custom_curve":
        exprs = params.get("exprs", {})
        t_sym = sp.symbols("t")
        try:
            fx = sp.lambdify(t_sym, sp.sympify(exprs.get("x", "0")), "numpy")
            fy = sp.lambdify(t_sym, sp.sympify(exprs.get("y", "0")), "numpy")
            fz = sp.lambdify(t_sym, sp.sympify(exprs.get("z", "0")), "numpy")
            X = np.array(fx(t), dtype=float).reshape(-1)
            Y = np.array(fy(t), dtype=float).reshape(-1)
            Z = np.array(fz(t), dtype=float).reshape(-1)
            R = np.column_stack((X, Y, Z))
        except Exception as e:
            raise ValueError(f"Error parsing custom curve expressions: {e}")

    else:
        raise ValueError(f"Unknown curve: {curve_name}")

    return t, R


# Surfaces


def get_surface_expressions(surface, params):
    if surface == "sphere":
        r = params.get("r", 1.0)

        return {
            "x": f"{r}*sin(u)*cos(v)",
            "y": f"{r}*sin(u)*sin(v)",
            "z": f"{r}*cos(u)",
        }

    if surface == "torus":
        R = params.get("R", 1.0)
        r = params.get("r", 0.4)
        return {
            "x": f"({R} + {r}*cos(u))*cos(v)",
            "y": f"({R} + {r}*cos(u))*sin(v)",
            "z": f"{r}*sin(u)",
        }

    if surface == "paraboloid":
        a = params.get("a", 1.0)
        return {
            "x": "u*cos(v)",
            "y": "u*sin(v)",
            "z": f"{a}*u**2",
        }

    # custom_surface
    return {
        "x": params.get("x", "0"),
        "y": params.get("y", "0"),
        "z": params.get("z", "0"),
    }


def mesh_from_parametric_surfaces(exprs, u_range, v_range, nu, nv):
    u = np.linspace(u_range[0], u_range[1], nu)
    v = np.linspace(v_range[0], v_range[1], nv)
    U, V = np.meshgrid(u, v, indexing="xy")

    usym, vsym = sp.symbols("u v")
    try:
        fx = sp.lambdify((usym, vsym), sp.sympify(exprs["x"]), "numpy")
        fy = sp.lambdify((usym, vsym), sp.sympify(exprs["y"]), "numpy")
        fz = sp.lambdify((usym, vsym), sp.sympify(exprs["z"]), "numpy")

        X = np.array(fx(U, V), dtype=float)
        Y = np.array(fy(U, V), dtype=float)
        Z = np.array(fz(U, V), dtype=float)

    except Exception as e:
        raise ValueError(f"Error evaluating surface expressions: {e}")

    return U, V, X, Y, Z


def symbolic_formula_for(curve, params):
    """
    Return a SymPy Matrix parametrization r(t) for the chosen curve.
    (Views.py will later convert this to LaTeX with sp.latex.)
    """
    t = sp.symbols("t")

    if curve == "line":
        P = sp.Matrix(
            [
                params.get("x0", 0),
                params.get("y0", 0),
                params.get("z0", 0),
            ]
        )
        Q = sp.Matrix(
            [
                params.get("x1", 1),
                params.get("y1", 0),
                params.get("z1", 0),
            ]
        )
        expr = P + t * (Q - P)

    elif curve == "circle":
        a = sp.symbols("a")
        expr = sp.Matrix([a * sp.cos(t), a * sp.sin(t)])

    elif curve == "ellipse":
        a, b = sp.symbols("a b")
        expr = sp.Matrix([a * sp.cos(t), b * sp.sin(t)])

    elif curve == "helix":
        a, b = sp.symbols("a b")
        expr = sp.Matrix([a * sp.cos(t), a * sp.sin(t), b * t])

    elif curve == "cycloid":
        a = sp.symbols("a")
        expr = sp.Matrix([a * (t - sp.sin(t)), a * (1 - sp.cos(t))])

    elif curve == "twisted_cubic":
        expr = sp.Matrix([t, t**2, t**3])

    elif curve == "catenary":
        C = sp.symbols("C")
        expr = sp.Matrix([t, C * sp.cosh(t / C)])

    elif curve == "hyperbola":
        expr = sp.Matrix([sp.cosh(t), sp.sinh(t)])

    elif curve == "tractrix":
        expr = sp.Matrix([t - sp.tanh(t), sp.sech(t)])

    else:
        raise ValueError(f"Unknown curve type: {curve}")

    return expr


def find_constants(parametrization, parameters):
    """
    Identify constants in a given parametrization
    """
    if not isinstance(parameters, list):
        parameters = [parameters]

    # Collect all symbols that are variables
    param_set = set(parameters)

    quantity_set = set()
    symbol_set = set()
    constants = []
    for expr in parametrization:
        if not isinstance(expr, Expr):
            continue  # skip ints, floats, etc.

        quantity_set |= expr.atoms(Quantity)
        symbol_set |= expr.free_symbols

    for symbol in symbol_set:
        if str(symbol) not in [str(q) for q in param_set]:
            constant = sp.Symbol(str(symbol), real=True, positive=True)
            constants.append(constant)

    return constants


def trig_simplify_all(expr):
    """
    Attempt a bunch of trig simplifications and return
    the simplest one because sp is a mongoloid apparently
    """
    if not isinstance(expr, sp.Expr):
        expr = sp.sympify(expr)

    candidates = set()

    # Basic simplifications
    candidates.add(sp.simplify(expr))
    candidates.add(sp.trigsimp(expr))
    candidates.add(sp.powsimp(expr))

    # Expansions
    candidates.add(sp.expand(expr))
    candidates.add(sp.expand_trig(expr))
    candidates.add(sp.expand_complex(expr))

    # Factorizations
    candidates.add(sp.factor(expr))
    candidates.add(sp.factor(expr.rewrite(sp.sin)))
    candidates.add(sp.factor(expr.rewrite(sp.cos)))

    # Rewrite in various trig bases
    for fn in (sp.sin, sp.cos, sp.tan, sp.exp):
        candidates.add(sp.trigsimp(expr.rewrite(fn)))
        candidates.add(sp.simplify(expr.rewrite(fn)))

    # Eliminate common subexpressions
    cse_expr, _ = sp.cse(expr)
    if cse_expr:
        candidates.add(cse_expr[-1])

    # Remove None and non-Expr
    candidates = [c for c in candidates if isinstance(c, sp.Expr)]

    # Choose the simplest by sp size metric
    best = min(candidates, key=lambda e: e.count_ops())

    return best


def true_simplify(expr):
    """
    You would think something like this is unnecessary in 2025 but alas
    """
    expr = sp.simplify(expr)
    syms = [s for s in expr.free_symbols]
    new_syms = {s: sp.Symbol(str(s), real=True, positive=True) for s in syms}

    expr_new = expr.subs(new_syms)
    expr_new = sp.simplify(expr_new)
    # expr_new = trig_simplify_all(expr_new)
    expr_new = expr_new.replace(lambda x: isinstance(x, sp.Abs), lambda x: x.args[0])

    return sp.simplify(expr_new)


def compute_unit_normal_vector(parametrization, parameters):
    """
    Calculate the unit normal vector of a surface given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    Matrix
        The unit normal vector as a sp Matrix.
    """

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])

    # Step 2 - Compute the normal vector using cross product
    normal_vector = X_u.cross(X_v)

    # Step 3 - Normalize the normal vector to get the unit normal vector
    magnitude = normal_vector.norm()

    magnitude = sp.simplify(magnitude)
    magnitude_no_abs = magnitude.replace(sp.Abs, lambda x: x)

    unit_normal = normal_vector / magnitude_no_abs

    # Step 4 - Simplify
    unit_normal = sp.simplify(unit_normal)

    return unit_normal


def compute_arc_length(parametrization, parameter, bounds):
    """
    Calculate the arc length of a curve given its parametrization

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the curve
    parameter : sp symbol
        A sp symbol representing the parameter of the curve
    bounds : list
        A list specifying the integration bounds for the parameter
    Returns:
    sp expression
        The arc length ds
    """
    parametrization = sp.Matrix(parametrization)

    # Step 1 - Compute the derivative of the parametrization
    X_t = sp.Matrix([parametrization.diff(parameter)])

    # Step 2 - Compute the magnitude of the derivative
    magnitude = X_t.norm()

    # Step 3 - Simplify
    magnitude = sp.simplify(magnitude)
    magnitude = magnitude.replace(sp.Abs, lambda x: x)

    # Step 4 - Integrate to get the arc length
    arc_length_integral = sp.Integral(magnitude, (parameter, bounds[0], bounds[1]))

    try:  # In case the integral is too hard
        arc_length = func_timeout(6, arc_length_integral.doit)
    except FunctionTimedOut:
        dict = {
            "msg": "No elementary antiderivative found.",
            "integral": arc_length_integral,
        }
        return dict

    # Step 5 - Simplify again
    arc_length = sp.simplify(arc_length)

    return arc_length


def compute_arc_length_reparametrization(parametrization, parameter, bounds):
    """
    Reparametrize a curve by arc length

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the curve.
    parameter : sp symbol
        A sp symbol representing the parameter of the curve.
    bounds : list
        A list specifying the integration bounds for the parameter.
    Returns:
    sp expression
        The arc length element ds.
    """
    # Step 1 - Get arc length
    arc_length = compute_arc_length(parametrization, parameter, bounds)

    if isinstance(arc_length, dict):
        return arc_length  # Return the dict with message if integral failed

    # Step 2 - Solve for t in terms of s
    s_dummy = sp.Dummy("s", real=True)
    equation = sp.Eq(arc_length, s_dummy)

    t_in_terms_of_s = sp.solve(equation, parameter)[0]

    # Step 3 - Substitute t back into the parametrization
    reparametrized_curve = [
        coord.subs(parameter, t_in_terms_of_s) for coord in parametrization
    ]

    # Step 4 - Simplify
    reparametrized_curve = [sp.simplify(coord) for coord in reparametrized_curve]

    return reparametrized_curve


def compute_first_fundamental_form(parametrization, parameters):
    """
    Calculate the first fundamental form of a surface given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
        Matrix
        The first fundamental form matrix.
    """

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])

    # Step 2 - Compute the coefficients of the first fundamental form using dot product
    E = X_u.dot(X_u)
    F = X_u.dot(X_v)
    G = X_v.dot(X_v)

    # Step 3 - Simplify
    E = sp.simplify(E)
    F = sp.simplify(F)
    G = sp.simplify(G)

    # Step 4 - Construct the first fundamental form matrix
    first_fundamental_form_matrix = Matrix([[E, F], [F, G]])

    return first_fundamental_form_matrix


def compute_surface_area(parametrization, parameters, u_bounds, v_bounds):
    """
    Calculate the surface area element of a surface given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.
    u_bounds : list
        A list specifying the integration bounds for the u parameter.
    v_bounds : list
        A list specifying the integration bounds for the v parameter.
    Returns:
    sp expression
        The surface area element dA.
    """

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )

    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute the integrand
    integrand = sp.refine(
        sp.sqrt(E * G - F**2), Q.positive(E) & Q.positive(G) & Q.positive(E * G - F**2)
    )
    integrand_no_abs = integrand.replace(sp.Abs, lambda x: x)

    # Step 3 - Compute the surface area
    surface_area_integral = sp.Integral(
        integrand_no_abs,
        (parameters[0], u_bounds[0], u_bounds[1]),
        (parameters[1], v_bounds[0], v_bounds[1]),
    )

    # In case the integral is too hard
    try:
        surface_area = func_timeout(6, surface_area_integral.doit)
    except FunctionTimedOut:
        dict = {
            "msg": "No elementary antiderivative found.",
            "integral": surface_area_integral,
        }
        return dict

    return surface_area


def compute_frenet_serret_apparatus(parametrization, parameter):
    """
    Calculate the Frenet-Serret apparatus of a curve given its parametrization.

    Parameters:
    parametrization : list
        A list of sympy expressions representing the parametrization
        of the curve, e.g. [x(t), y(t)] or [x(t), y(t), z(t)].
    parameter : sympy symbol
        The parameter of the curve (usually t).

    Returns:
    dict
        {
          "T": unit tangent vector,
          "N": unit normal vector,
          "B": binormal vector,
          "kappa": curvature,
          "tau": torsion
        }
    """

    # Step 1 - First derivative (velocity) and simplify
    X_t = Matrix([sp.diff(coord, parameter) for coord in parametrization])
    X_t = sp.simplify(X_t)

    # Step 2 - Unit tangent vector T
    X_t_norm = sp.sqrt(X_t.dot(X_t))
    T = X_t / X_t_norm
    T = sp.simplify(T)

    # Step 3 - Second derivative (acceleration)
    X_tt = Matrix([sp.diff(coord, parameter) for coord in X_t])
    X_tt = sp.simplify(X_tt)

    if all(sp.simplify(coord) == 0 for coord in X_tt):
        raise ValueError("Frenet-Serret apparatus undefined for straight lines.")

    # Step 4 - Curvature kappa
    kappa = (X_t.cross(X_tt)).norm() / (X_t_norm**3)
    kappa = sp.simplify(kappa)

    # Step 5 - Normal vector N
    T_t = Matrix([sp.diff(comp, parameter) for comp in T])
    T_t_norm = sp.sqrt(T_t.dot(T_t))
    N = sp.simplify(T_t / T_t_norm)

    # Step 6 - Binormal B
    B = sp.simplify(T.cross(N))

    # Step 7 - Torsion tau
    X_ttt = Matrix([sp.diff(coord, parameter) for coord in X_tt])
    tau = (X_t.cross(X_tt)).dot(X_ttt) / (X_t.cross(X_tt)).norm() ** 2
    tau = sp.simplify(tau)

    frenet_serret_dict = {"T": T, "N": N, "B": B, "kappa": kappa, "tau": tau}
    return frenet_serret_dict


def compute_second_fundamental_form(parametrization, parameters):
    """
    Calculate the second fundamental form of a surface given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    Matrix
        The second fundamental form matrix.
    """

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])
    X_uu = Matrix([sp.diff(coord, parameters[0]) for coord in X_u])
    X_uv = Matrix([sp.diff(coord, parameters[1]) for coord in X_u])
    X_vv = Matrix([sp.diff(coord, parameters[1]) for coord in X_v])

    # Step 2 - Compute the unit normal vector
    N = compute_unit_normal_vector(parametrization, parameters)

    # Step 3 - Compute the coefficients of the second fundamental form using dot product
    L = N.dot(X_uu)
    M = N.dot(X_uv)
    N_coeff = N.dot(X_vv)

    # Step 4 - Simplify
    L = sp.simplify(L)
    M = sp.simplify(M)
    N_coeff = sp.simplify(N_coeff)

    # Step 5 - Construct the second fundamental form matrix
    second_fundamental_form_matrix = Matrix([[L, M], [M, N_coeff]])

    return second_fundamental_form_matrix


def compute_shape_operator(parametrization, parameters):
    """
    Calculate the shape operator of a surface given its parametrization

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    Matrix
        The shape operator matrix.
    """
    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )

    # Step 2 - Compute the second fundamental form matrix
    second_fundamental_form_matrix = compute_second_fundamental_form(
        parametrization, parameters
    )

    # Step 3 - Compute the inverse of the first fundamental form matrix
    if first_fundamental_form_matrix.det() == 0:
        raise ValueError("Shape operator undefined for degenerate parametrization.")
    I_inv = first_fundamental_form_matrix.inv()
    I_inv = sp.simplify(I_inv)

    # Step 4 - Compute the shape operator as the product of I_inv and II
    shape_operator_matrix = I_inv * second_fundamental_form_matrix

    # Step 5 - Simplify
    shape_operator_matrix = sp.simplify(shape_operator_matrix)

    return shape_operator_matrix


def compute_mean_curvature(parametrization, parameters):
    """
    Calculate the mean curvature of a surface given its parametrization

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    sp expression
        The mean curvature H.
    """
    # Step 1 - Compute the shape operator matrix
    shape_operator_matrix = compute_shape_operator(parametrization, parameters)

    # Step 2 - Compute mean curvature as half the trace of the shape operator
    H = (shape_operator_matrix[0, 0] + shape_operator_matrix[1, 1]) / 2

    # Step 3 - Simplify
    H = sp.simplify(H)

    return H


def compute_gaussian_curvature(parametrization, parameters):
    """
    Calculate the Gaussian curvature of a surface given its parametrization

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    sp expression
        The Gaussian curvature K.
    """
    # Step 1 - Compute the shape operator matrix
    shape_operator_matrix = compute_shape_operator(parametrization, parameters)

    # Step 2 - Compute Gaussian curvature using determinant of shape operator
    K = shape_operator_matrix.det()

    # Step 3 - Simplify
    K = sp.simplify(K)

    return K


def compute_christoffel_symbols(parametrization, parameters):
    """
    Calculate the Christoffel symbols of the first kind from the first fundamental form.

    Parameters:
    first_fundamental_form : Matrix
        The first fundamental form matrix.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    dict
        A dictionary with keys as tuples (i, j, k) representing the Christoffel symbols Γ^k_ij.
    """

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )
    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute partial derivatives
    E_u = sp.diff(E, parameters[0])
    E_v = sp.diff(E, parameters[1])
    F_u = sp.diff(F, parameters[0])
    F_v = sp.diff(F, parameters[1])
    G_u = sp.diff(G, parameters[0])
    G_v = sp.diff(G, parameters[1])

    # Step 3 - Compute the inverse of the first fundamental form matrix and multipliers for each pair of Christoffel symbols
    I_inv = first_fundamental_form_matrix.inv()

    half = sp.Rational(1, 2)

    gamma_sub_uu_multiplier = Matrix([[half * E_u], [F_u - half * E_v]])
    gamma_sub_uv_multiplier = Matrix([[half * E_v], [half * G_u]])
    gamma_sub_vv_multiplier = Matrix([[F_v - half * G_u], [half * G_v]])

    # Step 4 - Compute the Christoffel symbols with matrix multiplication
    gamma_u_sub_uu = I_inv.row(0) * gamma_sub_uu_multiplier
    gamma_v_sub_uu = I_inv.row(1) * gamma_sub_uu_multiplier
    gamma_u_sub_uv = I_inv.row(0) * gamma_sub_uv_multiplier
    gamma_v_sub_uv = I_inv.row(1) * gamma_sub_uv_multiplier
    gamma_u_sub_vv = I_inv.row(0) * gamma_sub_vv_multiplier
    gamma_v_sub_vv = I_inv.row(1) * gamma_sub_vv_multiplier

    gamma = {
        "gamma_u_sub_uu": gamma_u_sub_uu,
        "gamma_v_sub_uu": gamma_v_sub_uu,
        "gamma_u_sub_uv": gamma_u_sub_uv,
        "gamma_v_sub_uv": gamma_v_sub_uv,
        "gamma_u_sub_vv": gamma_u_sub_vv,
        "gamma_v_sub_vv": gamma_v_sub_vv,
    }

    # Step 5 - Simplify
    for key, expr in gamma.items():
        gamma[key] = sp.simplify(expr[0])

    return gamma


def compute_codazzi_equations(parametrization, parameters):
    """
    Calculate the Codazzi equations of the surface.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    tuple
        A tuple containing the two Codazzi equations.
    """
    # Step 1 - Compute the second fundamental form matrix
    second_fundamental_form_matrix = compute_second_fundamental_form(
        parametrization, parameters
    )
    L = second_fundamental_form_matrix[0, 0]
    M = second_fundamental_form_matrix[0, 1]
    N_coeff = second_fundamental_form_matrix[1, 1]

    # Step 2 - Compute the Christoffel symbols
    christoffel_symbols = compute_christoffel_symbols(parametrization, parameters)
    gamma_u_sub_uu = dict(christoffel_symbols)["gamma_u_sub_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["gamma_v_sub_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["gamma_u_sub_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["gamma_v_sub_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["gamma_u_sub_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["gamma_v_sub_vv"]

    # Step 3 - Compute the Codazzi equations
    first_codazzi_eq_rhs = (
        L * gamma_u_sub_uv
        + M * gamma_v_sub_uv
        - M * gamma_u_sub_uu
        - N_coeff * gamma_v_sub_uu
    )
    second_codazzi_eq_rhs = (
        L * gamma_u_sub_vv
        + M * gamma_v_sub_vv
        - M * gamma_u_sub_uv
        - N_coeff * gamma_v_sub_uv
    )

    # Step 4 - Simplify
    first_codazzi_eq_rhs = sp.simplify(first_codazzi_eq_rhs)
    second_codazzi_eq_rhs = sp.simplify(second_codazzi_eq_rhs)

    # Step 5 - Verify
    L_v = sp.diff(L, parameters[1])
    M_u = sp.diff(M, parameters[0])
    M_v = sp.diff(M, parameters[1])
    N_u = sp.diff(N_coeff, parameters[0])

    first_codazzi_eq_lhs = sp.simplify(L_v - M_u)
    second_codazzi_eq_lhs = sp.simplify(M_v - N_u)

    return (
        first_codazzi_eq_rhs,
        second_codazzi_eq_rhs,
        first_codazzi_eq_lhs,
        second_codazzi_eq_lhs,
    )


def has_float(expr):
    return any(isinstance(a, sp.Float) for a in expr.atoms(sp.Float))


def compute_gauss_equations(parametrization, parameters):
    """
    Calculate the Gauss equations of the surface.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.

    Returns:
    tuple
        A tuple containing the four Gauss equations.
    """
    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )
    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute the second fundamental form matrix
    second_fundamental_form_matrix = compute_second_fundamental_form(
        parametrization, parameters
    )
    L = second_fundamental_form_matrix[0, 0]
    M = second_fundamental_form_matrix[0, 1]
    N = second_fundamental_form_matrix[1, 1]

    # Step 3 - Compute the Christoffel symbols
    christoffel_symbols = compute_christoffel_symbols(parametrization, parameters)
    gamma_u_sub_uu = dict(christoffel_symbols)["gamma_u_sub_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["gamma_v_sub_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["gamma_u_sub_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["gamma_v_sub_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["gamma_u_sub_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["gamma_v_sub_vv"]

    # Step 4 - Compute partial derivatives of Christoffel symbols
    gamma_u_sub_uv_u = sp.diff(gamma_u_sub_uv, parameters[0])
    gamma_u_sub_uu_v = sp.diff(gamma_u_sub_uu, parameters[1])
    gamma_v_sub_uv_u = sp.diff(gamma_v_sub_uv, parameters[0])
    gamma_v_sub_uu_v = sp.diff(gamma_v_sub_uu, parameters[1])
    gamma_u_sub_vv_u = sp.diff(gamma_u_sub_vv, parameters[0])
    gamma_u_sub_uv_v = sp.diff(gamma_u_sub_uv, parameters[1])
    gamma_v_sub_vv_u = sp.diff(gamma_v_sub_vv, parameters[0])
    gamma_v_sub_uv_v = sp.diff(gamma_v_sub_uv, parameters[1])

    # Step 5 - Compute the Gauss equations
    first_gauss_eq = (
        gamma_v_sub_uu_v
        - gamma_v_sub_uv_u
        + gamma_u_sub_uu * gamma_v_sub_uv
        + gamma_v_sub_uu * gamma_v_sub_vv
        - gamma_u_sub_uv * gamma_v_sub_uu
        - gamma_v_sub_uv**2
    )
    second_gauss_eq = (
        gamma_u_sub_uv_u
        - gamma_u_sub_uu_v
        + gamma_v_sub_uv * gamma_u_sub_uv
        - gamma_v_sub_uu * gamma_u_sub_vv
    )

    third_gauss_eq = (
        gamma_v_sub_uv_v
        - gamma_v_sub_vv_u
        + gamma_u_sub_uv * gamma_v_sub_uv
        - gamma_u_sub_vv * gamma_v_sub_uu
    )

    fourth_gauss_eq = (
        gamma_u_sub_vv_u
        - gamma_u_sub_uv_v
        + gamma_u_sub_vv * gamma_u_sub_uu
        - gamma_v_sub_vv * gamma_u_sub_uv
        - gamma_u_sub_uv**2
        - gamma_v_sub_uv * gamma_u_sub_vv
    )

    # Step 6 - Simplify

    first_gauss_eq_rhs = sp.simplify(first_gauss_eq)
    second_gauss_eq_rhs = sp.simplify(second_gauss_eq)
    third_gauss_eq_rhs = sp.simplify(third_gauss_eq)
    fourth_gauss_eq_rhs = sp.simplify(fourth_gauss_eq)

    K_numerator = L * N - M**2
    K_denominator = E * G - F**2
    K = sp.simplify(K_numerator / K_denominator)

    first_gauss_eq_lhs = sp.simplify(E * K)
    second_gauss_eq_lhs = sp.simplify(F * K)
    third_gauss_eq_lhs = sp.simplify(F * K)
    fourth_gauss_eq_lhs = sp.simplify(G * K)

    return (
        first_gauss_eq_rhs,
        second_gauss_eq_rhs,
        third_gauss_eq_rhs,
        fourth_gauss_eq_rhs,
        first_gauss_eq_lhs,
        second_gauss_eq_lhs,
        third_gauss_eq_lhs,
        fourth_gauss_eq_lhs,
    )
