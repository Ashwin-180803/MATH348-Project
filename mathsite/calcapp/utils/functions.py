import numpy as np
import sympy as sp
from sympy import Matrix, Q, integrate, sympify, Expr, assuming
from sympy.physics.units.quantities import Quantity
from func_timeout import func_timeout, FunctionTimedOut
import math
import os
import re


def parse_range_value(expr_str, default=0.0):
    if not expr_str or not isinstance(expr_str, str):
        return float(default)

    expr_str = expr_str.strip()
    if not expr_str:
        return float(default)

    try:
        return float(expr_str)
    except ValueError:
        try:
            expr = sp.sympify(expr_str)

            result = float(expr.evalf())
            return result
        except (ValueError, TypeError, AttributeError) as e:
            print(
                f"Warning: Could not parse range expression '{expr_str}', using default {default}: {e}"
            )
            return float(default)


# Curves


def get_default_curve_expressions(curve, params):
    """Get the default expressions for a curve type as strings, matching param_script.js defaults"""
    if curve == "helix":
        return {
            "x": "R * cos(t)",
            "y": "R * sin(t)",
            "z": "a * t",
        }

    if curve == "circle":
        return {
            "x": "R * cos(t)",
            "y": "R * sin(t)",
            "z": "0",
        }

    if curve == "ellipse":
        return {
            "x": "a * cos(t)",
            "y": "b * sin(t)",
            "z": "0",
        }

    if curve == "line":
        return {
            "x": "x0 + t*(x1-x0)",
            "y": "y0 + t*(y1-y0)",
            "z": "z0 + t*(z1-z0)",
        }

    if curve == "cycloid":
        return {
            "x": "r*(t - sin(t))",
            "y": "r*(1 - cos(t))",
            "z": "0",
        }

    if curve == "twisted_cubic":
        return {
            "x": "t",
            "y": "t**2",
            "z": "t**3",
        }

    if curve == "catenary":
        return {
            "x": "t",
            "y": "a * cosh(t/a)",
            "z": "0",
        }

    if curve == "hyperbola":
        return {
            "x": "cosh(t)",
            "y": "sinh(t)",
            "z": "0",
        }

    if curve == "tractrix":
        return {
            "x": "t - tanh(t)",
            "y": "sech(t)",
            "z": "0",
        }

    # Default fallback for custom_curve or unknown
    return {
        "x": "0",
        "y": "0",
        "z": "0",
    }


def numeric_curve_positions(curve_name, params, t0, t1, n):
    t = np.linspace(t0, t1, n)

    exprs = params.get("exprs", {})
    x_expr_str = str(exprs.get("x", "")).strip() if exprs.get("x") else ""
    y_expr_str = str(exprs.get("y", "")).strip() if exprs.get("y") else ""
    if exprs and x_expr_str and y_expr_str:
        t_sym = sp.symbols("t")
        try:
            z_expr_str = str(exprs.get("z", "0")).strip() if exprs.get("z") else "0"
            x_expr = sp.sympify(x_expr_str) if x_expr_str else sp.sympify("0")
            y_expr = sp.sympify(y_expr_str) if y_expr_str else sp.sympify("0")
            z_expr = sp.sympify(z_expr_str) if z_expr_str else sp.sympify("0")

            subs_dict = {}
            for key, value in params.items():
                if key != "exprs":
                    try:
                        param_value = float(value)
                        subs_dict[sp.Symbol(key)] = param_value
                    except (ValueError, TypeError):
                        pass

            if subs_dict:
                x_expr = x_expr.subs(subs_dict)
                y_expr = y_expr.subs(subs_dict)
                z_expr = z_expr.subs(subs_dict)

            fx = sp.lambdify(t_sym, x_expr, "numpy")
            fy = sp.lambdify(t_sym, y_expr, "numpy")
            fz = sp.lambdify(t_sym, z_expr, "numpy")

            X_raw = fx(t)
            Y_raw = fy(t)
            Z_raw = fz(t)

            X = np.asarray(X_raw, dtype=float).flatten()
            Y = np.asarray(Y_raw, dtype=float).flatten()
            Z = np.asarray(Z_raw, dtype=float).flatten()

            expected_len = len(t)

            if X.ndim == 0 or len(X) == 1:
                X = np.full(expected_len, float(X.flat[0]))
            else:
                X = (
                    X[:expected_len]
                    if len(X) >= expected_len
                    else np.pad(
                        X,
                        (0, expected_len - len(X)),
                        mode="constant",
                        constant_values=0,
                    )
                )

            if Y.ndim == 0 or len(Y) == 1:
                Y = np.full(expected_len, float(Y.flat[0]))
            else:
                Y = (
                    Y[:expected_len]
                    if len(Y) >= expected_len
                    else np.pad(
                        Y,
                        (0, expected_len - len(Y)),
                        mode="constant",
                        constant_values=0,
                    )
                )

            if Z.ndim == 0 or len(Z) == 1:
                Z = np.full(expected_len, float(Z.flat[0]))
            else:
                Z = (
                    Z[:expected_len]
                    if len(Z) >= expected_len
                    else np.pad(
                        Z,
                        (0, expected_len - len(Z)),
                        mode="constant",
                        constant_values=0,
                    )
                )

            min_len = min(len(X), len(Y), len(Z), len(t))
            X = X[:min_len]
            Y = Y[:min_len]
            Z = Z[:min_len]
            t = t[:min_len]

            R = np.column_stack((X, Y, Z))
            return t, R
        except Exception as e:
            raise ValueError(f"Error parsing curve expressions: {e}")

    if curve_name == "line":
        P = np.array(
            [
                float(params.get("x0", 0.0)),
                float(params.get("y0", 0.0)),
                float(params.get("z0", 0.0)),
            ],
            dtype=float,
        )
        Q = np.array(
            [
                float(params.get("x1", 1.0)),
                float(params.get("y1", 0.0)),
                float(params.get("z1", 0.0)),
            ],
            dtype=float,
        )
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


def get_default_surface_expressions(surface, params):
    """Get the default expressions for a surface type, matching param_script.js defaults"""
    if surface == "plane":
        return {
            "x": "u",
            "y": "v",
            "z": "a * u + b * v",
        }

    if surface == "cylinder":
        return {
            "x": "cos(u)",
            "y": "sin(u)",
            "z": "v",
        }

    if surface == "cone":
        return {
            "x": "u * cos(v)",
            "y": "u * sin(v)",
            "z": "u",
        }

    if surface == "paraboloid":
        return {
            "x": "u",
            "y": "v",
            "z": "u^2 + v^2",
        }

    if surface == "hyperbolic_paraboloid":
        return {
            "x": "u",
            "y": "v",
            "z": "u^2 - v^2",
        }

    if surface == "sphere":
        return {
            "x": "R * sin(u) * cos(v)",
            "y": "R * sin(u) * sin(v)",
            "z": "R * cos(u)",
        }

    if surface == "torus":
        return {
            "x": "(R + r * cos(v)) * cos(u)",
            "y": "(R + r * cos(v)) * sin(u)",
            "z": "r * sin(v)",
        }

    if surface == "helicoid":
        return {
            "x": "v * cos(u)",
            "y": "v * sin(u)",
            "z": "u",
        }

    if surface == "catenoid":
        return {
            "x": "a * cosh(v/a) * cos(u)",
            "y": "a * cosh(v/a) * sin(u)",
            "z": "v",
        }

    if surface == "mobius":
        return {
            "x": "(1 + v * cos(u/2)) * cos(u)",
            "y": "(1 + v * cos(u/2)) * sin(u)",
            "z": "v * sin(u/2)",
        }

    if surface == "klein":
        return {
            "x": "(cos(u) * (cos(u/2) * (sqrt(2)+cos(v)) + sin(u/2) * sin(v)))",
            "y": "(sin(u) * (cos(u/2) * (sqrt(2)+cos(v)) + sin(u/2) * sin(v)))",
            "z": "(sin(u/2) * (sqrt(2)+cos(v)) - cos(u/2) * sin(v))",
        }

    if surface == "enneper":
        return {
            "x": "u - (u^3)/3 + u*v^2",
            "y": "v - (v^3)/3 + v*u^2",
            "z": "u^2 - v^2",
        }

    # Default fallback
    return {
        "x": "u",
        "y": "v",
        "z": "0",
    }


def substitute_params_in_expr(expr_str, params, surface):
    if not expr_str:
        return expr_str

    substitutions = {}

    if surface == "torus":
        R_val = float(params.get("R", 1.0))
        r_val = float(params.get("r", 0.4))

        if re.search(r"\bR\b", expr_str):
            substitutions["R"] = str(R_val)
        if re.search(r"\br\b", expr_str):
            substitutions["r"] = str(r_val)
    elif surface == "sphere":
        r_val = float(params.get("r", 1.0))
        if re.search(r"\br\b", expr_str):
            substitutions["r"] = str(r_val)
    elif surface == "paraboloid":
        a_val = float(params.get("a", 1.0))
        if re.search(r"\ba\b", expr_str):
            substitutions["a"] = str(a_val)

    result = expr_str
    for param_name, param_value in substitutions.items():
        pattern = r"\b" + re.escape(param_name) + r"\b"
        result = re.sub(pattern, param_value, result)

    return result


def get_surface_expressions(surface, params):
    x_expr = params.get("x", "").strip()
    y_expr = params.get("y", "").strip()
    z_expr = params.get("z", "").strip()

    default_exprs = get_default_surface_expressions(surface, params)

    if (
        (x_expr and x_expr != "0" and x_expr != default_exprs["x"])
        or (y_expr and y_expr != "0" and y_expr != default_exprs["y"])
        or (z_expr and z_expr != "0" and z_expr != default_exprs["z"])
    ):
        x_expr = substitute_params_in_expr(x_expr, params, surface)
        y_expr = substitute_params_in_expr(y_expr, params, surface)
        z_expr = substitute_params_in_expr(z_expr, params, surface)
        return {
            "x": x_expr,
            "y": y_expr,
            "z": z_expr,
        }

    default_exprs_substituted = {
        "x": substitute_params_in_expr(default_exprs["x"], params, surface),
        "y": substitute_params_in_expr(default_exprs["y"], params, surface),
        "z": substitute_params_in_expr(default_exprs["z"], params, surface),
    }
    return default_exprs_substituted


def mesh_from_parametric_surfaces(
    exprs, u_range, v_range, nu, nv, var_u="u", var_v="v"
):
    u = np.linspace(u_range[0], u_range[1], nu)
    v = np.linspace(v_range[0], v_range[1], nv)
    U, V = np.meshgrid(u, v, indexing="xy")

    usym, vsym = sp.symbols(var_u + " " + var_v)
    try:
        print(
            f"DEBUG: Evaluating expressions with variables {var_u}, {var_v}: x='{exprs['x']}', y='{exprs['y']}', z='{exprs['z']}'"
        )
        fx = sp.lambdify((usym, vsym), sp.sympify(exprs["x"]), "numpy")
        fy = sp.lambdify((usym, vsym), sp.sympify(exprs["y"]), "numpy")
        fz = sp.lambdify((usym, vsym), sp.sympify(exprs["z"]), "numpy")

        X = np.array(fx(U, V), dtype=float)
        Y = np.array(fy(U, V), dtype=float)
        Z = np.array(fz(U, V), dtype=float)
        print(
            f"DEBUG: Arrays created: X.shape={X.shape}, Y.shape={Y.shape}, Z.shape={Z.shape}"
        )

    except Exception as e:
        print(f"DEBUG: Error in mesh generation: {e}")
        raise ValueError(f"Error evaluating surface expressions: {e}")

    return U, V, X, Y, Z


def symbolic_formula_for(curve, params):
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

    elif curve == "custom_curve":
        expr = sp.Matrix([0, 0])

    else:
        raise ValueError(f"Unknown curve type: {curve}")

    return expr


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

    steps = {}

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])
    steps["1"] = [X_u, X_v]

    # Step 2 - Compute the normal vector using cross product
    normal_vector = X_u.cross(X_v)
    steps["2"] = [normal_vector]

    # Step 3 - Normalize the normal vector to get the unit normal vector
    magnitude = normal_vector.norm()

    magnitude = sp.simplify(magnitude)
    magnitude_no_abs = magnitude.replace(sp.Abs, lambda x: x)

    unit_normal = normal_vector / magnitude_no_abs
    steps["3"] = [unit_normal, magnitude_no_abs]

    # Step 4 - Simplify
    unit_normal = sp.simplify(unit_normal)
    steps["4"] = [unit_normal]

    return unit_normal, steps


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
    steps = {}

    # Step 1 - Compute the derivative of the parametrization
    X_t = sp.Matrix([parametrization.diff(parameter)])
    steps["1"] = [X_t]
    
    # Step 2 - Compute the magnitude of the derivative and simplify
    magnitude = X_t.norm()

    magnitude = sp.simplify(magnitude)
    magnitude = magnitude.replace(sp.Abs, lambda x: x)
    steps["2"] = [magnitude]

    # Step 3 - Integrate magnitude to get the arc length and simplify
    arc_length_integral = sp.Integral(magnitude, (parameter, bounds[0], bounds[1]))
    
    try:  # In case the integral is too hard
        arc_length = func_timeout(6, arc_length_integral.doit)
    except FunctionTimedOut:
        dict = {
            "msg": "No elementary antiderivative found.",
            "integral": arc_length_integral,
        }
        return dict

    arc_length = sp.simplify(arc_length)
    steps["3"] = [arc_length_integral, arc_length]

    return arc_length, steps


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
    steps = {}
    bounds[1] = parameter

    # Step 1 - Compute arc length
    arc_length, _ = compute_arc_length(parametrization, parameter, bounds)
    
    if isinstance(arc_length, dict):
        return arc_length  # Return the dict with message if integral failed

    steps["1"] = [arc_length]

    # Step 2 - Solve for t in terms of s
    s_dummy = sp.Dummy("s", real=True)
    equation = sp.Eq(arc_length, s_dummy)
    t_in_terms_of_s = sp.solve(equation, parameter)[0]

    steps["2"] = [t_in_terms_of_s]

    # Step 3 - Substitute t back into the parametrization and simplify
    reparametrized_curve = [
        coord.subs(parameter, t_in_terms_of_s) for coord in parametrization
    ]
    reparametrized_curve = [sp.simplify(coord) for coord in reparametrized_curve]
    steps["3"] = [reparametrized_curve]

    return reparametrized_curve, steps


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
    steps = {}

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])
    steps["1"] = [X_u, X_v]

    # Step 2 - Compute the coefficients of the first fundamental form using dot product and simplify
    E = X_u.dot(X_u)
    F = X_u.dot(X_v)
    G = X_v.dot(X_v)

    E = sp.simplify(E)
    F = sp.simplify(F)
    G = sp.simplify(G)

    steps["2"] = [E, F, G]

    # Step 3 - Construct the first fundamental form matrix
    first_fundamental_form_matrix = Matrix([[E, F], [F, G]])
    steps["3"] = [first_fundamental_form_matrix]

    return first_fundamental_form_matrix, steps


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
    steps = {}

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix, _ = compute_first_fundamental_form(
        parametrization, parameters
    )
    steps["1"] = [first_fundamental_form_matrix]

    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    # Step 2 - Compute the integrand
    integrand = sp.refine(
        sp.sqrt(E * G - F**2), Q.positive(E) & Q.positive(G) & Q.positive(E * G - F**2)
    )
    integrand_no_abs = integrand.replace(sp.Abs, lambda x: x)
    steps["2"] = [integrand_no_abs]

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

    surface_area = sp.simplify(surface_area)
    steps["3"] = [surface_area_integral, surface_area]

    return surface_area, steps


def compute_frenet_serret_apparatus(parametrization, parameter, allow_straight_lines=False):
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
    steps = {}
    straight_line_message = "Frenet-Serret apparatus undefined for straight lines."

    # Step 1 - Find first derivative (velocity) and simplify
    X_t = Matrix([sp.diff(coord, parameter) for coord in parametrization])
    X_t = sp.simplify(X_t)
    steps["1"] = [X_t]

    # Step 2 - Find unit tangent vector T
    X_t_norm = sp.sqrt(X_t.dot(X_t))
    T = X_t / X_t_norm
    T = sp.simplify(T)
    steps["2"] = [T]

    # Step 3 - Find second derivative (acceleration)
    X_tt = Matrix([sp.diff(coord, parameter) for coord in X_t])
    X_tt = sp.simplify(X_tt)
    

    if all(sp.simplify(coord) == 0 for coord in X_tt):
        if not allow_straight_lines:
            raise ValueError(straight_line_message)

        # Straight line fallback: curvature is zero, normal/binormal are undefined
        zero_curvature = sp.Integer(0)
        nan_vector = Matrix([sp.nan] * len(X_t))
        frenet_serret_dict = {
            "T": T,
            "N": nan_vector,
            "B": nan_vector,
            "kappa": zero_curvature,
            "tau": sp.nan,
        }
        steps["warning"] = [straight_line_message]
        return frenet_serret_dict, steps

    steps["3"] = [X_tt]

    # Step 4 - Find curvature kappa
    kappa = (X_t.cross(X_tt)).norm() / (X_t_norm**3)
    kappa = sp.simplify(kappa)
    steps["4"] = [kappa]

    # Step 5 - Find normal vector N
    T_t = Matrix([sp.diff(comp, parameter) for comp in T])
    T_t_norm = sp.sqrt(T_t.dot(T_t))
    N = sp.simplify(T_t / T_t_norm)
    steps["5"] = [N]

    # Step 6 - Find binormal B
    B = sp.simplify(T.cross(N))
    steps["6"] = [B]

    # Step 7 - Find torsion tau
    X_ttt = Matrix([sp.diff(coord, parameter) for coord in X_tt])
    tau = (X_t.cross(X_tt)).dot(X_ttt) / (X_t.cross(X_tt)).norm() ** 2
    tau = sp.simplify(tau)
    steps["7"] = [tau]

    frenet_serret_dict = {"T": T, "N": N, "B": B, "kappa": kappa, "tau": tau}
    return frenet_serret_dict, steps


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
    steps = {}

    # Step 1 - Compute the partial derivatives
    X_u = Matrix([sp.diff(coord, parameters[0]) for coord in parametrization])
    X_v = Matrix([sp.diff(coord, parameters[1]) for coord in parametrization])
    X_uu = Matrix([sp.diff(coord, parameters[0]) for coord in X_u])
    X_uv = Matrix([sp.diff(coord, parameters[1]) for coord in X_u])
    X_vv = Matrix([sp.diff(coord, parameters[1]) for coord in X_v])
    steps["1"] = [X_u, X_v, X_uu, X_uv, X_vv]

    # Step 2 - Compute the unit normal vector
    N, _ = compute_unit_normal_vector(parametrization, parameters)
    steps["2"] = [N]

    # Step 3 - Compute the coefficients of the second fundamental form using dot product
    L = N.dot(X_uu)
    M = N.dot(X_uv)
    N_coeff = N.dot(X_vv)

    L = sp.simplify(L)
    M = sp.simplify(M)
    N_coeff = sp.simplify(N_coeff)

    steps["3"] = [L, M, N_coeff]

    # Step 4 - Construct the second fundamental form matrix
    second_fundamental_form_matrix = Matrix([[L, M], [M, N_coeff]])
    steps["4"] = [second_fundamental_form_matrix]

    return second_fundamental_form_matrix, steps


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
    steps = {}

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix, _ = compute_first_fundamental_form(
        parametrization, parameters
    )
    steps["1"] = [first_fundamental_form_matrix]

    # Step 2 - Compute the second fundamental form matrix
    second_fundamental_form_matrix, _ = compute_second_fundamental_form(
        parametrization, parameters
    )
    steps["2"] = [second_fundamental_form_matrix]

    # Step 3 - Compute the inverse of the first fundamental form matrix
    if first_fundamental_form_matrix.det() == 0:
        raise ValueError("Shape operator undefined for degenerate parametrization.")
    I_inv = first_fundamental_form_matrix.inv()
    I_inv = sp.simplify(I_inv)
    steps["3"] = [I_inv]

    # Step 4 - Compute the shape operator as the product of I_inv and II and simplify
    shape_operator_matrix = I_inv * second_fundamental_form_matrix
    shape_operator_matrix = sp.simplify(shape_operator_matrix)

    steps["4"] = [shape_operator_matrix]

    return shape_operator_matrix, steps


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
    steps = {}

    # Step 1 - Compute the shape operator matrix
    shape_operator_matrix, _ = compute_shape_operator(parametrization, parameters)
    steps["1"] = [shape_operator_matrix]

    # Step 2 - Compute mean curvature as half the trace of the shape operator and simplify
    H = (shape_operator_matrix[0, 0] + shape_operator_matrix[1, 1]) / 2
    H = sp.simplify(H)

    steps["2"] = [H]

    return H, steps


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
    steps = {}

    # Step 1 - Compute the shape operator matrix
    shape_operator_matrix, _ = compute_shape_operator(parametrization, parameters)
    steps["1"] = [shape_operator_matrix]

    # Step 2 - Compute Gaussian curvature using determinant of shape operator and simplify
    K = shape_operator_matrix.det()
    K = sp.simplify(K)
    steps["2"] = [K]

    return K, steps


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
    steps = {}

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix, _ = compute_first_fundamental_form(
        parametrization, parameters
    )
    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]
    steps["1"] = [first_fundamental_form_matrix]

    # Step 2 - Compute partial derivatives
    E_u = sp.diff(E, parameters[0])
    E_v = sp.diff(E, parameters[1])
    F_u = sp.diff(F, parameters[0])
    F_v = sp.diff(F, parameters[1])
    G_u = sp.diff(G, parameters[0])
    G_v = sp.diff(G, parameters[1])
    steps["2"] = [E_u, E_v, F_u, F_v, G_u, G_v]

    # Step 3 - Compute the inverse of the first fundamental form matrix and multipliers for each pair of Christoffel symbols
    I_inv = first_fundamental_form_matrix.inv()

    half = sp.Rational(1, 2)

    gamma_sub_uu_multiplier = Matrix([[half * E_u], [F_u - half * E_v]])
    gamma_sub_uv_multiplier = Matrix([[half * E_v], [half * G_u]])
    gamma_sub_vv_multiplier = Matrix([[F_v - half * G_u], [half * G_v]])

    steps["3"] = [
        I_inv,
        gamma_sub_uu_multiplier,
        gamma_sub_uv_multiplier,
        gamma_sub_vv_multiplier,
    ]

    # Step 4 - Compute the Christoffel symbols with matrix multiplication and simplify
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

    for key, expr in gamma.items():
        gamma[key] = sp.simplify(expr[0])

    steps["4"] = gamma

    return gamma, steps


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
    steps = {}

    # Step 1 - Compute the second fundamental form matrix
    second_fundamental_form_matrix, _ = compute_second_fundamental_form(
        parametrization, parameters
    )
    L = second_fundamental_form_matrix[0, 0]
    M = second_fundamental_form_matrix[0, 1]
    N_coeff = second_fundamental_form_matrix[1, 1]

    steps["1"] = [second_fundamental_form_matrix]

    # Step 2 - Compute the Christoffel symbols
    christoffel_symbols, _ = compute_christoffel_symbols(parametrization, parameters)
    gamma_u_sub_uu = dict(christoffel_symbols)["gamma_u_sub_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["gamma_v_sub_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["gamma_u_sub_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["gamma_v_sub_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["gamma_u_sub_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["gamma_v_sub_vv"]

    steps["2"] = christoffel_symbols

    # Step 3 - Compute the Codazzi equations and simplify
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

    first_codazzi_eq_rhs = sp.simplify(first_codazzi_eq_rhs)
    second_codazzi_eq_rhs = sp.simplify(second_codazzi_eq_rhs)

    steps["3"] = [first_codazzi_eq_rhs, second_codazzi_eq_rhs]

    # Step 4 - Verify
    L_v = sp.diff(L, parameters[1])
    M_u = sp.diff(M, parameters[0])
    M_v = sp.diff(M, parameters[1])
    N_u = sp.diff(N_coeff, parameters[0])

    first_codazzi_eq_lhs = sp.simplify(L_v - M_u)
    second_codazzi_eq_lhs = sp.simplify(M_v - N_u)

    steps["4"] = [first_codazzi_eq_lhs, second_codazzi_eq_lhs, L_v, M_u, M_v, N_u]

    return (
        first_codazzi_eq_rhs,
        second_codazzi_eq_rhs,
        first_codazzi_eq_lhs,
        second_codazzi_eq_lhs,
        steps,
    )


def has_float(expr):
    return any(isinstance(a, sp.Float) for a in expr.atoms(sp.Float))


def compute_numeric_frenet_serret(t, R):
    """
    Compute curvature and torsion numerically from position data.

    Parameters:
    t : array
        Parameter values
    R : array
        Position vectors [x, y, z] for each t

    Returns:
    dict
        {
            "curvature": array of curvature values,
            "torsion": array of torsion values,
            "arc_length": array of arc length values
        }
    """
    n = len(t)
    curvature = np.full(n, None)
    torsion = np.full(n, None)
    arc_length = np.zeros(n)

    # Compute arc length
    for i in range(1, n):
        dist = np.linalg.norm(R[i] - R[i - 1])
        arc_length[i] = arc_length[i - 1] + dist

    # Compute derivatives numerically using finite differences
    dt = np.gradient(t)

    # First derivatives (velocity)
    dR_dt = np.gradient(R, t, axis=0)

    # Speed (magnitude of velocity)
    speed = np.linalg.norm(dR_dt, axis=1)

    # Avoid division by zero
    speed = np.where(speed < 1e-10, 1e-10, speed)

    # Unit tangent vector T
    T = dR_dt / speed[:, np.newaxis]

    # Second derivatives (acceleration)
    dT_dt = np.gradient(T, t, axis=0)

    curvature_magnitude = np.linalg.norm(dT_dt, axis=1)
    curvature = curvature_magnitude / speed

    return {
        "curvature": curvature.tolist(),
        "torsion": torsion.tolist(),
        "arc_length": arc_length.tolist(),
    }


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
    steps = {}

    # Step 1 - Compute the first fundamental form matrix
    first_fundamental_form_matrix, _ = compute_first_fundamental_form(
        parametrization, parameters
    )
    E = first_fundamental_form_matrix[0, 0]
    F = first_fundamental_form_matrix[0, 1]
    G = first_fundamental_form_matrix[1, 1]

    steps["1"] = [first_fundamental_form_matrix]

    # Step 2 - Compute the second fundamental form matrix
    second_fundamental_form_matrix, _ = compute_second_fundamental_form(
        parametrization, parameters
    )
    L = second_fundamental_form_matrix[0, 0]
    M = second_fundamental_form_matrix[0, 1]
    N = second_fundamental_form_matrix[1, 1]

    steps["2"] = [second_fundamental_form_matrix]

    # Step 3 - Compute the Christoffel symbols
    christoffel_symbols, _ = compute_christoffel_symbols(parametrization, parameters)
    gamma_u_sub_uu = dict(christoffel_symbols)["gamma_u_sub_uu"]
    gamma_v_sub_uu = dict(christoffel_symbols)["gamma_v_sub_uu"]
    gamma_u_sub_uv = dict(christoffel_symbols)["gamma_u_sub_uv"]
    gamma_v_sub_uv = dict(christoffel_symbols)["gamma_v_sub_uv"]
    gamma_u_sub_vv = dict(christoffel_symbols)["gamma_u_sub_vv"]
    gamma_v_sub_vv = dict(christoffel_symbols)["gamma_v_sub_vv"]

    steps["3"] = christoffel_symbols

    # Step 4 - Compute partial derivatives of Christoffel symbols
    gamma_u_sub_uv_u = sp.diff(gamma_u_sub_uv, parameters[0])
    gamma_u_sub_uu_v = sp.diff(gamma_u_sub_uu, parameters[1])
    gamma_v_sub_uv_u = sp.diff(gamma_v_sub_uv, parameters[0])
    gamma_v_sub_uu_v = sp.diff(gamma_v_sub_uu, parameters[1])
    gamma_u_sub_vv_u = sp.diff(gamma_u_sub_vv, parameters[0])
    gamma_u_sub_uv_v = sp.diff(gamma_u_sub_uv, parameters[1])
    gamma_v_sub_vv_u = sp.diff(gamma_v_sub_vv, parameters[0])
    gamma_v_sub_uv_v = sp.diff(gamma_v_sub_uv, parameters[1])

    steps["4"] = [
        gamma_u_sub_uv_u,
        gamma_u_sub_uu_v,
        gamma_v_sub_uv_u,
        gamma_v_sub_uu_v,
        gamma_u_sub_vv_u,
        gamma_u_sub_uv_v,
        gamma_v_sub_vv_u,
        gamma_v_sub_uv_v,
    ]

    # Step 5 - Compute the Gauss equations and simplify
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

    first_gauss_eq_rhs = sp.simplify(first_gauss_eq)
    second_gauss_eq_rhs = sp.simplify(second_gauss_eq)
    third_gauss_eq_rhs = sp.simplify(third_gauss_eq)
    fourth_gauss_eq_rhs = sp.simplify(fourth_gauss_eq)

    steps["5"] = [
        first_gauss_eq_rhs,
        second_gauss_eq_rhs,
        third_gauss_eq_rhs,
        fourth_gauss_eq_rhs,
    ]

    # Step 6 - Compute Gaussian curvature K to verify

    K_numerator = L * N - M**2
    K_denominator = E * G - F**2
    K = sp.simplify(K_numerator / K_denominator)

    first_gauss_eq_lhs = sp.simplify(E * K)
    second_gauss_eq_lhs = sp.simplify(F * K)
    third_gauss_eq_lhs = sp.simplify(F * K)
    fourth_gauss_eq_lhs = sp.simplify(G * K)

    steps["6"] = [
        K,
        first_gauss_eq_lhs,
        second_gauss_eq_lhs,
        third_gauss_eq_lhs,
        fourth_gauss_eq_lhs,
    ]

    return (
        first_gauss_eq_rhs,
        second_gauss_eq_rhs,
        third_gauss_eq_rhs,
        fourth_gauss_eq_rhs,
        first_gauss_eq_lhs,
        second_gauss_eq_lhs,
        third_gauss_eq_lhs,
        fourth_gauss_eq_lhs,
        steps,
    )


def get_computation_steps(computation_type, parametrization=None, parameters=None):
    steps = []

    if computation_type == "arc_length":
        steps = [
            "Compute the derivative of the parametrization",
            "Compute the magnitude of the derivative",
            "Integrate magnitude to get the arc length",
        ]

    elif computation_type == "reparam_arc_length":
        steps = [
            "Compute arc length",
            "Solve for t in terms of s",
            "Substitute t back into the parametrization",
        ]

    elif computation_type == "frenet":
        steps = [
            "Find first derivative (velocity)",
            "Find unit tangent vector T",
            "Find second derivative (acceleration)",
            "Find curvature kappa",
            "Find normal vector N",
            "Find binormal B",
            "Find torsion tau",
        ]

    elif computation_type == "first_form":
        steps = [
            "Compute the partial derivatives",
            "Compute the coefficients of the first fundamental form using dot product",
            "Construct the first fundamental form matrix",
        ]

    elif computation_type == "second_form":
        steps = [
            "Compute the partial derivatives",
            "Compute the unit normal vector",
            "Compute the coefficients of the second fundamental form using dot product",
            "Construct the second fundamental form matrix",
        ]

    elif computation_type == "gaussian_curvature":
        steps = [
            "Compute the shape operator matrix",
            "Compute Gaussian curvature using determinant of shape operator",
        ]

    elif computation_type == "mean_curvature":
        steps = [
            "Compute the shape operator matrix",
            "Compute mean curvature as half the trace of the shape operator",
        ]

    elif computation_type == "principal_curvatures":
        steps = [
            "Compute the first fundamental form matrix",
            "Compute the second fundamental form matrix",
            "Compute the inverse of the first fundamental form matrix",
            "Compute the shape operator as the product of I_inv and II",
        ]

    elif computation_type == "christoffel":
        steps = [
            "Compute the first fundamental form matrix",
            "Compute partial derivatives",
            "Compute the inverse of the first fundamental form matrix and multipliers for each pair of Christoffel symbols",
            "Compute the Christoffel symbols with matrix multiplication",
        ]

    elif computation_type == "gauss_equations":
        steps = [
            "Compute the first fundamental form matrix",
            "Compute the second fundamental form matrix",
            "Compute the Christoffel symbols",
            "Compute partial derivatives of Christoffel symbols",
            "Compute the Gauss equations",
        ]

    elif computation_type == "codazzi_equations":
        steps = [
            "Compute the second fundamental form matrix",
            "Compute the Christoffel symbols",
            "Compute the Codazzi equations",
        ]

    return steps
