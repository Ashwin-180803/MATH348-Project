from flask import Flask, request, jsonify, send_from_directory
import numpy as np
import sympy as sp
from sympy import Matrix, Q
from sympy.physics.units.quantities import Quantity
import math
import os

app = Flask(__name__, static_folder="static", static_url_path="/static")


def numeric_curve_positions(curve_name, params, t0, t1, n):
    t = np.linspace(t0, t1, n)
    if curve_name == "line":
        P = np.array([params["x0"], params["y0"], params.get("z0", 0.0)], dtype=float)
        Q = np.array([params["x1"], params["y1"], params.get("z1", 0.0)], dtype=float)
        R = P[np.newaxis, :] + np.outer(t, (Q - P))
    elif curve_name == "circle":
        a = float(params.get("a", 1.0))
        R = np.column_stack((a * np.cos(t), a * np.sin(t), np.zeros_like(t)))
    elif curve_name == "ellipse":
        a = float(params.get("a", 2.0))
        b = float(params.get("b", 1.0))
        R = np.column_stack((a * np.cos(t), b * np.sin(t), np.zeros_like(t)))
    elif curve_name == "helix":
        a = float(params.get("a", 1.0))
        b = float(params.get("b", 0.2))
        R = np.column_stack((a * np.cos(t), a * np.sin(t), b * t))
    elif curve_name == "cycloid":
        a = float(params.get("a", 1.0))
        R = np.column_stack(
            (a * (t - np.sin(t)), a * (1 - np.cos(t)), np.zeros_like(t))
        )
    elif curve_name == "twisted_cubic":
        R = np.column_stack((t, t**2, t**3))
    elif curve_name == "catenary":
        C = float(params.get("C", 1.0))
        R = np.column_stack((t, C * np.cosh(t / C), np.zeros_like(t)))
    elif curve_name == "hyperbola":
        R = np.column_stack((np.cosh(t), np.sinh(t), np.zeros_like(t)))
    elif curve_name == "tractrix":
        R = np.column_stack((t - np.tanh(t), 1 / np.cosh(t), np.zeros_like(t)))
    else:
        raise ValueError("Unknown curve: " + curve_name)
    return t, R


def derivatives_numeric(R, t):
    dt = np.gradient(t)
    Rp = np.gradient(R, axis=0) / dt[:, np.newaxis]
    Rpp = np.gradient(Rp, axis=0) / dt[:, np.newaxis]
    Rppp = np.gradient(Rpp, axis=0) / dt[:, np.newaxis]
    return Rp, Rpp, Rppp


def curvature_torsion_from_R(R, t):
    Rp, Rpp, Rppp = derivatives_numeric(R, t)
    cross = np.cross(Rp, Rpp)
    cross_norm = np.linalg.norm(cross, axis=1)
    Rp_norm = np.linalg.norm(Rp, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        curvature = cross_norm / (Rp_norm**3)
    triple = np.einsum("ij,ij->i", cross, Rppp)
    denom = cross_norm**2
    with np.errstate(divide="ignore", invalid="ignore"):
        torsion = triple / denom
    curvature = [None if (not np.isfinite(x)) else float(x) for x in curvature]
    torsion = [None if (not np.isfinite(x)) else float(x) for x in torsion]
    return curvature, torsion


def arc_length_numeric(R, t):
    Rp = np.gradient(R, axis=0) / np.gradient(t)[:, np.newaxis]
    speed = np.linalg.norm(Rp, axis=1)
    dt = np.diff(t)
    if len(dt) == 0:
        return [0.0]
    cum = np.zeros_like(speed)
    cum[1:] = np.cumsum((speed[:-1] + speed[1:]) * 0.5 * dt)
    return cum.tolist()


@app.route("/")
def index():
    return send_from_directory("static", "index.html")


@app.route("/api/compute", methods=["POST"])
def api_compute():
    data = request.json
    curve = data.get("curve")
    params = data.get("params", {})
    t0 = float(data.get("t0", 0.0))
    t1 = float(data.get("t1", 2 * math.pi))
    n = int(data.get("n", 400))
    try:
        t, R = numeric_curve_positions(curve, params, t0, t1, n)
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 400

    curvature, torsion = curvature_torsion_from_R(R, t)
    arclen = arc_length_numeric(R, t)

    xs = R[:, 0].tolist()
    ys = R[:, 1].tolist()
    zs = R[:, 2].tolist()

    sym_formula = symbolic_formula_for(curve, params)

    return jsonify(
        {
            "ok": True,
            "curve": curve,
            "t": t.tolist(),
            "x": xs,
            "y": ys,
            "z": zs,
            "curvature": curvature,
            "torsion": torsion,
            "arc_length": arclen,
            "symbolic": sym_formula,
        }
    )


def symbolic_formula_for(curve, params):
    t = sp.symbols("t")
    if curve == "line":
        P = sp.Matrix([params.get("x0", 0), params.get("y0", 0), params.get("z0", 0)])
        Q = sp.Matrix([params.get("x1", 1), params.get("y1", 0), params.get("z1", 0)])
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
        return ""
    return sp.latex(expr)


def find_constants(parametrization, parameters):
    """
    Identify constants in a given parametrization
    """
    # Collect all symbols that are variables
    param_set = set(parameters)

    quantity_set = set().union(*(expr.atoms(Quantity) for expr in parametrization))

    symbol_set = set().union(*(expr.free_symbols for expr in parametrization))

    # True constants = quantities + non-variable symbols
    constants = quantity_set | (symbol_set - param_set)

    return constants


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
    u, v = parameters

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, u) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, v) for coord in parametrization])

    # Step 2 - Compute the normal vector using cross product
    normal_vector = X_u.cross(X_v)

    # Step 3 - Normalize the normal vector to get the unit normal vector
    magnitude = normal_vector.norm()
    unit_normal = normal_vector / magnitude

    # Step 4 - Simplify
    unit_normal = sp.simplify(unit_normal)

    return unit_normal


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


def compute_arc_length(parametrization, parameter, bounds):
    """
    Calculate the arc length element of a curve given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the curve.
    parameter : sp symbol
        A sp symbol representing the parameter of the curve.
    bounds : dict
        A dictionary specifying the integration bounds for the parameter.
    Returns:
    sp expression
        The arc length element ds.
    """
    t = parameter

    # Step 1 - Compute the derivative of the parametrization
    X_t = Matrix([sp.diff(coord, t) for coord in parametrization])

    # Step 2 - Compute the magnitude of the derivative
    magnitude = X_t.norm()

    # Step 3 - Simplify
    magnitude = trig_simplify_all(magnitude)

    # Step 4 - Integrate to get the arc length
    arc_length = sp.integrate(magnitude, (t, bounds[t][0], bounds[t][1]))

    return arc_length


def compute_arc_length_reparametrization(parametrization, parameter, bounds):
    """
    Reparametrize a curve by arc length

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the curve.
    parameter : sp symbol
        A sp symbol representing the parameter of the curve.
    bounds : dict
        A dictionary specifying the integration bounds for the parameter.
    Returns:
    sp expression
        The arc length element ds.
    """
    t = parameter

    # Step 1 - Get arc length
    arc_length = compute_arc_length(parametrization, parameter, bounds)

    # Step 2 - Solve for t in terms of s
    s = sp.symbols("s", real=True)
    equation = sp.Eq(arc_length, s)
    t_in_terms_of_s = sp.solve(equation, t)[0]

    # Step 3 - Substitute t back into the parametrization
    reparametrized_curve = [coord.subs(t, t_in_terms_of_s) for coord in parametrization]

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
    u, v = parameters

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, u) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, v) for coord in parametrization])

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


def compute_surface_area(parametrization, parameters, bounds):
    """
    Calculate the surface area element of a surface given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the surface.
    parameters : list
        A list of sp symbols representing the parameters of the surface.
    bounds : dict
        A dictionary specifying the integration bounds for each parameter.
    assumptions : dict, optional
        A dictionary of assumptions for the parameters, used to simplify the expressions.
    Returns:
    sp expression
        The surface area element dA.
    """
    u, v = parameters

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )

    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute the integrand
    # Make sure to get rid of absolute values if possible
    # I have to do it manually because sp is legitimately stupid
    with sp.assuming(Q.positive(E) & Q.positive(G) & Q.nonnegative(E * G - F**2)):
        integrand = sp.sqrt(E * G - F**2)
        integrand_no_abs = integrand.replace(
            lambda x: isinstance(x, sp.Abs), lambda x: x.args[0]
        )
        integrand_no_abs = sp.refine(
            integrand_no_abs,
            Q.positive(E) & Q.positive(G) & Q.nonnegative(E * G - F**2),
        )

    # Step 3 - Compute the surface area

    surface_area = sp.integrate(
        integrand_no_abs,
        (v, bounds[v][0], bounds[v][1]),
        (u, bounds[u][0], bounds[u][1]),
    )

    return surface_area


def compute_frenet_serret_apparatus(parametrization, parameter):
    """
    Calculate the Frenet-Serret apparatus of a curve given its parametrization.

    Parameters:
    parametrization : list
        A list of sp expressions representing the parametrization of the curve.
    parameter : sp symbol
        A sp symbol representing the parameter of the curve.

    Returns:
    tuple
        A tuple containing the tangent vector T, normal vector N, binormal vector B,
        curvature kappa, and torsion tau.
    """
    t = parameter

    # Step 1 - Compute the first derivative (velocity) and simplify
    X_t = Matrix([sp.diff(coord, t) for coord in parametrization])
    X_t = sp.simplify(X_t)

    # Step 2 - Compute the unit tangent vector T and simplify
    X_t_norm_factored = sp.factor(X_t.norm())
    T = X_t / X_t_norm_factored
    T = sp.simplify(T)

    # Step 3 - Compute the second derivative (acceleration) and simplify
    X_tt = Matrix([sp.diff(coord, t) for coord in X_t])
    X_tt = sp.simplify(X_tt)

    # Step 4 - Compute the curvature kappa and simplify
    kappa = (X_t.cross(X_tt)).norm() / (X_t.norm() ** 3)
    kappa = sp.simplify(kappa)

    # Step 5 - Compute the unit normal vector N and simplify
    N = sp.simplify(sp.diff(T, t) / sp.diff(T, t).norm())
    N = sp.simplify(N)

    # Step 6 - Compute the binormal vector B and simplify
    B = T.cross(N)
    B = sp.simplify(B)

    # Step 7 - Compute the torsion tau and simplify
    X_ttt = Matrix([sp.diff(coord, t) for coord in X_tt])
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
    u, v = parameters

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, u) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, v) for coord in parametrization])
    X_uu = Matrix([sp.diff(coord, u) for coord in X_u])
    X_uv = Matrix([sp.diff(coord, v) for coord in X_u])
    X_vv = Matrix([sp.diff(coord, v) for coord in X_v])

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
    I_inv = first_fundamental_form_matrix.inv()

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
    u, v = parameters

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix = compute_first_fundamental_form(
        parametrization, parameters
    )
    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute partial derivatives
    E_u = sp.diff(E, u)
    E_v = sp.diff(E, v)
    F_u = sp.diff(F, u)
    F_v = sp.diff(F, v)
    G_u = sp.diff(G, u)
    G_v = sp.diff(G, v)

    # Step 3 - Compute the inverse of the first fundamental form matrix and multipliers for each pair of Christoffel symbols
    I_inv = first_fundamental_form_matrix.inv()

    gamma_sub_uu_multiplier = Matrix([[0.5 * E_u], [F_u - 0.5 * E_v]])
    gamma_sub_uv_multiplier = Matrix([[0.5 * E_v], [0.5 * G_u]])
    gamma_sub_vv_multiplier = Matrix([[F_v - 0.5 * G_u], [0.5 * G_v]])

    # Step 4 - Compute the Christoffel symbols with matrix multiplication
    gamma_u_sub_uu = I_inv.row(0) * gamma_sub_uu_multiplier
    gamma_v_sub_uu = I_inv.row(1) * gamma_sub_uu_multiplier
    gamma_u_sub_uv = I_inv.row(0) * gamma_sub_uv_multiplier
    gamma_v_sub_uv = I_inv.row(1) * gamma_sub_uv_multiplier
    gamma_u_sub_vv = I_inv.row(0) * gamma_sub_vv_multiplier
    gamma_v_sub_vv = I_inv.row(1) * gamma_sub_vv_multiplier

    gamma = {
        "Γ^u_uu": gamma_u_sub_uu,
        "Γ^v_uu": gamma_v_sub_uu,
        "Γ^u_uv": gamma_u_sub_uv,
        "Γ^v_uv": gamma_v_sub_uv,
        "Γ^u_vv": gamma_u_sub_vv,
        "Γ^v_vv": gamma_v_sub_vv,
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
    gamma_u_sub_uu = dict(christoffel_symbols)["Γ^u_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["Γ^v_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["Γ^u_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["Γ^v_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["Γ^u_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["Γ^v_vv"]

    # Step 3 - Compute the Codazzi equations
    first_codazzi_eq = (
        L * gamma_u_sub_uv
        + M * gamma_v_sub_uv
        - M * gamma_u_sub_uu
        - N_coeff * gamma_v_sub_uu
    )
    second_codazzi_eq = (
        L * gamma_u_sub_vv
        + M * gamma_v_sub_vv
        - M * gamma_u_sub_uv
        - N_coeff * gamma_v_sub_uv
    )

    # Step 4 - Simplify
    first_codazzi_eq = sp.simplify(first_codazzi_eq)
    second_codazzi_eq = sp.simplify(second_codazzi_eq)

    L_v = sp.diff(L, parameters[1])
    M_u = sp.diff(M, parameters[0])
    M_v = sp.diff(M, parameters[1])
    N_u = sp.diff(N_coeff, parameters[0])

    return first_codazzi_eq, second_codazzi_eq


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
    gamma_u_sub_uu = dict(christoffel_symbols)["Γ^u_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["Γ^v_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["Γ^u_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["Γ^v_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["Γ^u_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["Γ^v_vv"]

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
    first_gauss_eq = sp.simplify(first_gauss_eq)
    second_gauss_eq = sp.simplify(second_gauss_eq)
    third_gauss_eq = sp.simplify(third_gauss_eq)
    fourth_gauss_eq = sp.simplify(fourth_gauss_eq)

    K_numerator = L * N - M**2
    K_denominator = E * G - F**2
    K = sp.simplify(K_numerator / K_denominator)

    return (
        first_gauss_eq,
        second_gauss_eq,
        third_gauss_eq,
        fourth_gauss_eq,
    )


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
