import numpy as np
import sympy as sp
from sympy import Matrix, Q, integrate, sympify, Expr, assuming
from sympy.physics.units.quantities import Quantity
from func_timeout import func_timeout, FunctionTimedOut
import math
import os
import re
from . import parser


def substitute_constants_in_expr(expr, params, parameters):
    
    if not params or not isinstance(expr, sp.Expr):
        return expr
    
    free_symbols = expr.free_symbols
    
    param_names = {str(p) for p in parameters}
    
    subs_dict = {}
    for key, value in params.items():
        if key in param_names:
            continue
        
        try:
            param_value = float(value)
        except (ValueError, TypeError):
            continue
        
        for sym in free_symbols:
            if str(sym) == key:
                subs_dict[sym] = param_value
                break
    
    if subs_dict:
        return expr.subs(subs_dict)
    return expr


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
<<<<<<< HEAD
            expr = parser.parse_input(expr_str, [])
            
=======
            expr = sp.sympify(expr_str)

>>>>>>> origin
            result = float(expr.evalf())
            return result
        except (ValueError, TypeError, AttributeError) as e:
            print(
                f"Warning: Could not parse range expression '{expr_str}', using default {default}: {e}"
            )
            return float(default)


# Curves


def get_default_curve_expressions(curve, params):
    
    if curve == "helix":
        return {
            "x": "cos(t)",
            "y": "sin(t)",
            "z": "0.2*t",
        }

    if curve == "circle":
        return {
            "x": "cos(t)",
            "y": "sin(t)",
            "z": "0",
        }

    if curve == "ellipse":
        return {
            "x": "2*cos(t)",
            "y": "sin(t)",
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
            "x": "t - sin(t)",
            "y": "1 - cos(t)",
            "z": "0",
        }

    if curve == "twisted_cubic":
        return {
            "x": "t",
            "y": "t^2",
            "z": "t^3",
        }

    if curve == "catenary":
        return {
            "x": "t",
            "y": "cosh(t/2)",
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
<<<<<<< HEAD
    
    
=======

    # Default fallback for custom_curve or unknown
>>>>>>> origin
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
<<<<<<< HEAD
            x_expr = parser.parse_input(x_expr_str, [t_sym]) if x_expr_str else parser.parse_input("0", [t_sym])
            y_expr = parser.parse_input(y_expr_str, [t_sym]) if y_expr_str else parser.parse_input("0", [t_sym])
            z_expr = parser.parse_input(z_expr_str, [t_sym]) if z_expr_str else parser.parse_input("0", [t_sym])
            
            x_expr = substitute_constants_in_expr(x_expr, params, [t_sym])
            y_expr = substitute_constants_in_expr(y_expr, params, [t_sym])
            z_expr = substitute_constants_in_expr(z_expr, params, [t_sym])
            
=======
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

>>>>>>> origin
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
        if not isinstance(exprs, dict):
            exprs = {}
        t_sym = sp.symbols("t")
        try:
            x_expr_str = str(exprs.get("x", "0")).strip() if exprs.get("x") else "0"
            y_expr_str = str(exprs.get("y", "0")).strip() if exprs.get("y") else "0"
            z_expr_str = str(exprs.get("z", "0")).strip() if exprs.get("z") else "0"
            
            if not x_expr_str:
                x_expr_str = "0"
            if not y_expr_str:
                y_expr_str = "0"
            if not z_expr_str:
                z_expr_str = "0"
            
            x_expr = parser.parse_input(x_expr_str, [t_sym])
            y_expr = parser.parse_input(y_expr_str, [t_sym])
            z_expr = parser.parse_input(z_expr_str, [t_sym])
            
            x_expr = substitute_constants_in_expr(x_expr, params, [t_sym])
            y_expr = substitute_constants_in_expr(y_expr, params, [t_sym])
            z_expr = substitute_constants_in_expr(z_expr, params, [t_sym])
            
            fx = sp.lambdify(t_sym, x_expr, "numpy")
            fy = sp.lambdify(t_sym, y_expr, "numpy")
            fz = sp.lambdify(t_sym, z_expr, "numpy")
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
            "z": "u + v",
        }

    if surface == "cylinder":
        return {
            "x": "cos(u)",
            "y": "sin(u)",
            "z": "v",
        }

    if surface == "cone":
        return {
            "x": "v * cos(u)",
            "y": "v * sin(u)",
            "z": "v",
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
            "x": "r * cos(u) * sin(v)",
            "y": "r * sin(u) * sin(v)",
            "z": "r * cos(v)",
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
            "x": "cosh(v) * cos(u)",
            "y": "cosh(v) * sin(u)",
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


def substitute_curve_params_in_expr(expr_str, params, curve):
    
    import re

    if not expr_str:
        return expr_str
<<<<<<< HEAD
    
    substitutions = {}
    
    
    if curve == "catenary":
        C_val = float(params.get("C", 1.0))
        if re.search(r'\bC\b', expr_str):
            substitutions["C"] = str(C_val)
    
    elif curve == "circle":
        a_val = float(params.get("a", 1.0))
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif curve == "ellipse":
        a_val = float(params.get("a", 2.0))
        b_val = float(params.get("b", 1.0))
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
        if re.search(r'\bb\b', expr_str):
            substitutions["b"] = str(b_val)
    
    elif curve == "helix":
        a_val = float(params.get("a", 1.0))
        b_val = float(params.get("b", 0.2))
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
        if re.search(r'\bb\b', expr_str):
            substitutions["b"] = str(b_val)
    
    elif curve == "cycloid":
        a_val = float(params.get("a", 1.0))
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif curve == "line":
        x0_val = float(params.get("x0", 0.0))
        y0_val = float(params.get("y0", 0.0))
        z0_val = float(params.get("z0", 0.0))
        x1_val = float(params.get("x1", 1.0))
        y1_val = float(params.get("y1", 0.0))
        z1_val = float(params.get("z1", 0.0))
        if re.search(r'\bx0\b', expr_str):
            substitutions["x0"] = str(x0_val)
        if re.search(r'\by0\b', expr_str):
            substitutions["y0"] = str(y0_val)
        if re.search(r'\bz0\b', expr_str):
            substitutions["z0"] = str(z0_val)
        if re.search(r'\bx1\b', expr_str):
            substitutions["x1"] = str(x1_val)
        if re.search(r'\by1\b', expr_str):
            substitutions["y1"] = str(y1_val)
        if re.search(r'\bz1\b', expr_str):
            substitutions["z1"] = str(z1_val)
    
    # Apply substitutions
    result = expr_str
    for param_name, param_value in substitutions.items():
        pattern = r'\b' + re.escape(param_name) + r'\b'
        result = re.sub(pattern, param_value, result)
    
    return result


def substitute_params_in_expr(expr_str, params, surface):
    
    import re
    if not expr_str:
        return expr_str
    
    substitutions = {}
    
    
    if surface == "torus":
        R_str = params.get("R", "").strip()
        r_str = params.get("r", "").strip()
        try:
            R_val = float(R_str) if R_str else 1.0
        except (ValueError, TypeError):
            R_val = 1.0
        try:
            r_val = float(r_str) if r_str else 0.4
        except (ValueError, TypeError):
            r_val = 0.4
        
        if re.search(r'\bR\b', expr_str):
=======

    substitutions = {}

    if surface == "torus":
        R_val = float(params.get("R", 1.0))
        r_val = float(params.get("r", 0.4))

        if re.search(r"\bR\b", expr_str):
>>>>>>> origin
            substitutions["R"] = str(R_val)
        if re.search(r"\br\b", expr_str):
            substitutions["r"] = str(r_val)
    
    elif surface == "sphere":
<<<<<<< HEAD
        r_str = params.get("r", "").strip()
        try:
            r_val = float(r_str) if r_str else 1.0
        except (ValueError, TypeError):
            r_val = 1.0
        if re.search(r'\br\b', expr_str):
=======
        r_val = float(params.get("r", 1.0))
        if re.search(r"\br\b", expr_str):
>>>>>>> origin
            substitutions["r"] = str(r_val)
    
    elif surface == "paraboloid":
<<<<<<< HEAD
        a_str = params.get("a", "").strip()
        try:
            a_val = float(a_str) if a_str else 1.0
        except (ValueError, TypeError):
            a_val = 1.0
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif surface == "hyperbolic_paraboloid":
        a_str = params.get("a", "").strip()
        try:
            a_val = float(a_str) if a_str else 1.0
        except (ValueError, TypeError):
            a_val = 1.0
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif surface == "cylinder":
        r_str = params.get("r", "").strip()
        try:
            r_val = float(r_str) if r_str else 1.0
        except (ValueError, TypeError):
            r_val = 1.0
        if re.search(r'\br\b', expr_str):
            substitutions["r"] = str(r_val)
    
    elif surface == "cone":
        r_str = params.get("r", "").strip()
        h_str = params.get("h", "").strip()
        try:
            r_val = float(r_str) if r_str else 1.0
        except (ValueError, TypeError):
            r_val = 1.0
        try:
            h_val = float(h_str) if h_str else 1.0
        except (ValueError, TypeError):
            h_val = 1.0
        if re.search(r'\br\b', expr_str):
            substitutions["r"] = str(r_val)
        if re.search(r'\bh\b', expr_str):
            substitutions["h"] = str(h_val)
    
    elif surface == "catenoid":
        a_str = params.get("a", "").strip()
        try:
            a_val = float(a_str) if a_str else 1.0
        except (ValueError, TypeError):
            a_val = 1.0
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif surface == "mobius":
        r_str = params.get("r", "").strip()
        try:
            r_val = float(r_str) if r_str else 1.0
        except (ValueError, TypeError):
            r_val = 1.0
        if re.search(r'\br\b', expr_str):
            substitutions["r"] = str(r_val)
    
    elif surface == "klein":
        
        a_str = params.get("a", "").strip()
        try:
            a_val = float(a_str) if a_str else 1.0
        except (ValueError, TypeError):
            a_val = 1.0
        if re.search(r'\ba\b', expr_str):
            substitutions["a"] = str(a_val)
    
    elif surface == "enneper":
       
        a_str = params.get("a", "").strip()
        try:
            a_val = float(a_str) if a_str else 1.0
        except (ValueError, TypeError):
            a_val = 1.0
        if re.search(r'\ba\b', expr_str):
=======
        a_val = float(params.get("a", 1.0))
        if re.search(r"\ba\b", expr_str):
>>>>>>> origin
            substitutions["a"] = str(a_val)

    result = expr_str
    for param_name, param_value in substitutions.items():
<<<<<<< HEAD
        pattern = r'\b' + re.escape(param_name) + r'\b'
=======
        pattern = r"\b" + re.escape(param_name) + r"\b"
>>>>>>> origin
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


<<<<<<< HEAD
def mesh_from_parametric_surfaces(exprs, u_range, v_range, nu, nv, var_u="u", var_v="v", params=None):
    """
    Generate mesh from parametric surface expressions.
    
    Args:
        exprs: Dictionary with "x", "y", "z" expression strings
        u_range: Tuple (u0, u1) for u parameter range
        v_range: Tuple (v0, v1) for v parameter range
        nu: Number of points in u direction
        nv: Number of points in v direction
        var_u: Variable name for u parameter (default "u")
        var_v: Variable name for v parameter (default "v")
        params: Optional dictionary of parameters for constant substitution (e.g., {"R": "1.0", "r": "0.4"})
    """
=======
def mesh_from_parametric_surfaces(
    exprs, u_range, v_range, nu, nv, var_u="u", var_v="v"
):
>>>>>>> origin
    u = np.linspace(u_range[0], u_range[1], nu)
    v = np.linspace(v_range[0], v_range[1], nv)
    U, V = np.meshgrid(u, v, indexing="xy")

    usym, vsym = sp.symbols(var_u + " " + var_v)
    try:
<<<<<<< HEAD
        print(f"DEBUG: Evaluating expressions with variables {var_u}, {var_v}: x='{exprs['x']}', y='{exprs['y']}', z='{exprs['z']}'")
        
        x_sympy = parser.parse_input(exprs["x"], [usym, vsym])
        y_sympy = parser.parse_input(exprs["y"], [usym, vsym])
        z_sympy = parser.parse_input(exprs["z"], [usym, vsym])
        
        if params:
            x_sympy = substitute_constants_in_expr(x_sympy, params, [usym, vsym])
            y_sympy = substitute_constants_in_expr(y_sympy, params, [usym, vsym])
            z_sympy = substitute_constants_in_expr(z_sympy, params, [usym, vsym])
        
        fx = sp.lambdify((usym, vsym), x_sympy, "numpy")
        fy = sp.lambdify((usym, vsym), y_sympy, "numpy")
        fz = sp.lambdify((usym, vsym), z_sympy, "numpy")

        
        X_raw = fx(U, V)
        Y_raw = fy(U, V)
        Z_raw = fz(U, V)
        
        
        X = np.asarray(X_raw, dtype=float)
        Y = np.asarray(Y_raw, dtype=float)
        Z = np.asarray(Z_raw, dtype=float)
        
        
        if X.ndim == 0:
            X = np.full_like(U, float(X))
        elif X.ndim == 1:
            X = np.tile(X, (len(v), 1))
        
        if Y.ndim == 0:
            Y = np.full_like(U, float(Y))
        elif Y.ndim == 1:
            Y = np.tile(Y, (len(v), 1))
        
        if Z.ndim == 0:
            Z = np.full_like(U, float(Z))
        elif Z.ndim == 1:
            Z = np.tile(Z, (len(v), 1))
        
       
        if X.shape != U.shape:
            X = np.broadcast_to(X, U.shape)
        if Y.shape != U.shape:
            Y = np.broadcast_to(Y, U.shape)
        if Z.shape != U.shape:
            Z = np.broadcast_to(Z, U.shape)
        
        print(f"DEBUG: Arrays created: X.shape={X.shape}, Y.shape={Y.shape}, Z.shape={Z.shape}")
        
        
        if X.size == 0 or Y.size == 0 or Z.size == 0:
            raise ValueError(f"Generated arrays are empty: X.size={X.size}, Y.size={Y.size}, Z.size={Z.size}")
=======
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
>>>>>>> origin

    except Exception as e:
        print(f"DEBUG: Error in mesh generation: {e}")
        import traceback
        print(f"DEBUG: Traceback: {traceback.format_exc()}")
        raise ValueError(f"Error evaluating surface expressions: {e}")

    return U, V, X, Y, Z


def symbolic_formula_for(curve, params):
<<<<<<< HEAD
    """
    Get symbolic formula for a curve, with parameter values substituted.
    Returns SymPy Matrix with numeric values substituted for parameters.
    """
=======
>>>>>>> origin
    t = sp.symbols("t")

    if curve == "line":
        x0_val = float(params.get("x0", 0.0))
        y0_val = float(params.get("y0", 0.0))
        z0_val = float(params.get("z0", 0.0))
        x1_val = float(params.get("x1", 1.0))
        y1_val = float(params.get("y1", 0.0))
        z1_val = float(params.get("z1", 0.0))
        P = sp.Matrix([x0_val, y0_val, z0_val])
        Q = sp.Matrix([x1_val, y1_val, z1_val])
        expr = P + t * (Q - P)

    elif curve == "circle":
        a_val = float(params.get("a", 1.0))
        expr = sp.Matrix([a_val * sp.cos(t), a_val * sp.sin(t)])

    elif curve == "ellipse":
        a_val = float(params.get("a", 2.0))
        b_val = float(params.get("b", 1.0))
        expr = sp.Matrix([a_val * sp.cos(t), b_val * sp.sin(t)])

    elif curve == "helix":
        a_val = float(params.get("a", 1.0))
        b_val = float(params.get("b", 0.2))
        expr = sp.Matrix([a_val * sp.cos(t), a_val * sp.sin(t), b_val * t])

    elif curve == "cycloid":
        a_val = float(params.get("a", 1.0))
        expr = sp.Matrix([a_val * (t - sp.sin(t)), a_val * (1 - sp.cos(t))])

    elif curve == "twisted_cubic":
        expr = sp.Matrix([t, t**2, t**3])

    elif curve == "catenary":
        C_val = float(params.get("C", 1.0))
        expr = sp.Matrix([t, C_val * sp.cosh(t / C_val)])

    elif curve == "hyperbola":
        expr = sp.Matrix([sp.cosh(t), sp.sinh(t)])

    elif curve == "tractrix":
        expr = sp.Matrix([t - sp.tanh(t), sp.sech(t)])

    else:
        raise ValueError(f"Unknown curve type: {curve}")

    return expr


def trig_simplify_all(expr):
    """
    Attempt a bunch of trig simplifications and return
    the simplest one because sp is a mongoloid apparently
    """
    if not isinstance(expr, sp.Expr):
        if isinstance(expr, str):
            try:
                expr = parser.parse_input(expr, [])
            except:
                expr = sp.sympify(expr)
        else:
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

    # Step 3 - Check if normal vector is zero (degenerate case)
    normal_vector_simplified = sp.simplify(normal_vector)
    if all(sp.simplify(comp) == 0 for comp in normal_vector_simplified):
        raise ValueError("Unit normal vector undefined: partial derivatives are linearly dependent (degenerate parametrization).")

    # Step 4 - Normalize the normal vector to get the unit normal vector
    magnitude = normal_vector.norm()
    magnitude = sp.simplify(magnitude)
    
    # Check if magnitude is zero
    magnitude_simplified = sp.simplify(magnitude)
    try:
        if magnitude_simplified == 0 or sp.simplify(magnitude_simplified - 0) == 0:
            raise ValueError("Unit normal vector undefined: cross product of partial derivatives is zero (degenerate parametrization).")
        if magnitude_simplified.is_zero:
            raise ValueError("Unit normal vector undefined: cross product of partial derivatives is zero (degenerate parametrization).")
    except (AttributeError, TypeError):
        # Can't determine symbolically, try to proceed
        pass
    
    magnitude_no_abs = magnitude.replace(sp.Abs, lambda x: x.args[0] if len(x.args) > 0 else x)

    try:
        unit_normal = normal_vector / magnitude_no_abs
    except (ZeroDivisionError, ValueError) as e:
        raise ValueError(f"Unit normal vector undefined: cannot normalize (magnitude is zero or near zero). This indicates a degenerate parametrization. {str(e)}")

    # Step 5 - Simplify
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
        
        # Check if the result is still an Integral (evaluation failed silently)
        if isinstance(arc_length, sp.Integral):
            # The integral couldn't be evaluated
            return {
                "msg": "No elementary antiderivative found.",
                "integral": arc_length_integral,
            }
    except FunctionTimedOut:
        return {
            "msg": "No elementary antiderivative found.",
            "integral": arc_length_integral,
        }
    except Exception as e:
        # If evaluation fails for other reasons, return error
        return {
            "msg": f"Error evaluating arc length integral: {str(e)}",
            "integral": arc_length_integral,
        }

    # Step 5 - Simplify again
    arc_length = sp.simplify(arc_length)
    
    # Final check: if it still contains an Integral, return error
    if isinstance(arc_length, sp.Integral) or arc_length.has(sp.Integral):
        return {
            "msg": "No elementary antiderivative found.",
            "integral": arc_length_integral,
        }

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
        # If arc length computation failed, we can't do reparametrization
        # Enhance the error message to be more informative
        if "integral" in arc_length:
            return {
                "msg": "Arc length reparametrization failed: The arc length integral for this curve does not have an elementary antiderivative. This means the arc length cannot be expressed in closed form using standard functions. Curves like ellipses, certain polynomials of degree 4 or higher, and curves with complex expressions often have this limitation.",
                "integral": arc_length.get("integral"),
                "error_type": "NoElementaryAntiderivative"
            }
        return arc_length  # Return the dict with message if integral failed

    # Step 1.5 - Check if arc_length contains an unevaluated Integral
    # This can happen if the integral evaluation timed out but returned the Integral object
    if isinstance(arc_length, sp.Integral) or arc_length.has(sp.Integral):
        # Try to evaluate the integral if it's still there
        try:
            # If it's a definite integral, try to evaluate it with a longer timeout
            arc_length_eval = func_timeout(10, lambda: arc_length.doit() if isinstance(arc_length, sp.Integral) else arc_length)
            if not isinstance(arc_length_eval, sp.Integral) and not arc_length_eval.has(sp.Integral):
                arc_length = arc_length_eval
                arc_length = sp.simplify(arc_length)
            else:
                # Still contains integral, can't proceed
                return {
                    "msg": "Cannot compute arc length reparametrization: the arc length integral could not be evaluated analytically. This curve may not have a closed-form arc length expression.",
                    "error_type": "IntegralNotEvaluated"
                }
        except FunctionTimedOut:
            return {
                "msg": "Cannot compute arc length reparametrization: the arc length integral evaluation timed out. This curve may not have a closed-form arc length expression.",
                "error_type": "IntegralTimeout"
            }
        except Exception as e:
            # If we can't evaluate it, we can't solve for reparametrization
            return {
                "msg": f"Cannot compute arc length reparametrization: the arc length integral could not be evaluated. {str(e)}",
                "error_type": "IntegralNotEvaluated"
            }
    
    # Check again if it still has Integral after evaluation attempt
    if isinstance(arc_length, sp.Integral) or arc_length.has(sp.Integral):
        return {
            "msg": "Cannot compute arc length reparametrization: the arc length involves an integral that cannot be evaluated analytically. This curve may not have a closed-form arc length expression.",
            "error_type": "IntegralNotEvaluated"
        }

    # Step 2 - Simplify arc length expression before solving
    arc_length_simplified = sp.simplify(arc_length)
    # Try to expand and simplify further
    try:
        arc_length_simplified = sp.expand(arc_length_simplified)
        arc_length_simplified = sp.simplify(arc_length_simplified)
    except:
        pass
    
    # Step 3 - Solve for t in terms of s
    try:
        s_dummy = sp.Dummy("s", real=True)
        s_sym = sp.symbols("s", real=True)
        
        # Try multiple solving methods
        solutions = None
        
        # Method 1: Try with simplified arc length
        try:
            # Check if we can create the equation (avoid Relational errors with Integrals)
            if isinstance(arc_length_simplified, sp.Integral) or arc_length_simplified.has(sp.Integral):
                raise ValueError("Arc length contains unevaluated integral")
            equation = sp.Eq(arc_length_simplified, s_sym)
            solutions = sp.solve(equation, parameter, dict=False)
            if solutions and len(solutions) > 0:
                # Filter out complex solutions if parameter is real
                real_solutions = [sol for sol in solutions if sol.is_real or not sol.has(sp.I)]
                if real_solutions:
                    solutions = real_solutions
        except Exception:
            pass
        
        # Method 2: Try with original arc length
        if not solutions or len(solutions) == 0:
            try:
                # Check if we can create the equation
                if isinstance(arc_length, sp.Integral) or arc_length.has(sp.Integral):
                    raise ValueError("Arc length contains unevaluated integral")
                equation = sp.Eq(arc_length, s_sym)
                solutions = sp.solve(equation, parameter, dict=False)
                if solutions and len(solutions) > 0:
                    real_solutions = [sol for sol in solutions if sol.is_real or not sol.has(sp.I)]
                    if real_solutions:
                        solutions = real_solutions
            except Exception:
                pass
        
        # Method 3: Try solve with manual=True (more thorough)
        if not solutions or len(solutions) == 0:
            try:
                # Check if we can create the equation
                if isinstance(arc_length_simplified, sp.Integral) or arc_length_simplified.has(sp.Integral):
                    raise ValueError("Arc length contains unevaluated integral")
                equation = sp.Eq(arc_length_simplified, s_sym)
                solutions = sp.solve(equation, parameter, dict=False, manual=True)
                if solutions and len(solutions) > 0:
                    real_solutions = [sol for sol in solutions if sol.is_real or not sol.has(sp.I)]
                    if real_solutions:
                        solutions = real_solutions
            except Exception:
                pass
        
        # Method 4: Try using the indefinite integral approach
        # Instead of solving s = ∫[t0 to t] ds/dt dt, try solving s = ∫ ds/dt dt (indefinite)
        if not solutions or len(solutions) == 0:
            try:
                # Compute ds/dt = ||X'(t)||
                parametrization_matrix = sp.Matrix(parametrization)
                X_t = sp.Matrix([parametrization_matrix.diff(parameter)])
                ds_dt = X_t.norm()
                ds_dt = sp.simplify(ds_dt)
                ds_dt = ds_dt.replace(sp.Abs, lambda x: x.args[0] if len(x.args) > 0 else x)
                
                # If speed is constant, t = s / speed + constant
                if ds_dt.is_constant() or ds_dt.is_number:
                    # If speed is constant, s = speed * (t - t0), so t = s/speed + t0
                    # For bounds starting at bounds[0], t = s/ds_dt + bounds[0]
                    t_solution = s_sym / ds_dt + bounds[0]
                    solutions = [t_solution]
                else:
                    # Try computing indefinite integral s(t) = ∫ ds/dt dt
                    try:
                        s_of_t = sp.integrate(ds_dt, parameter)
                        s_of_t = sp.simplify(s_of_t)
                        # Now solve s_of_t = s_sym for parameter
                        if not s_of_t.has(sp.Integral):
                            eq_indef = sp.Eq(s_of_t, s_sym)
                            solutions_indef = sp.solve(eq_indef, parameter, dict=False)
                            if solutions_indef and len(solutions_indef) > 0:
                                real_solutions = [sol for sol in solutions_indef if sol.is_real or not sol.has(sp.I)]
                                if real_solutions:
                                    solutions = real_solutions
                    except:
                        pass
            except Exception:
                pass
        
        # Method 5: Try rearranging and using different variable names
        if not solutions or len(solutions) == 0:
            try:
                # Check if we can create the equation
                if isinstance(arc_length_simplified, sp.Integral) or arc_length_simplified.has(sp.Integral):
                    raise ValueError("Arc length contains unevaluated integral")
                # Sometimes rewriting helps
                equation = sp.Eq(arc_length_simplified - s_sym, 0)
                solutions = sp.solve(equation, parameter, dict=False)
                if solutions and len(solutions) > 0:
                    real_solutions = [sol for sol in solutions if sol.is_real or not sol.has(sp.I)]
                    if real_solutions:
                        solutions = real_solutions
            except Exception:
                pass
        
        # Method 6: Try using nsolve for numerical approximation (as last resort)
        # But this requires a guess, so we'll skip it for now
        
        if not solutions or len(solutions) == 0:
            # Provide more helpful error message
            equation = None
            try:
                # Check if we can create the equation (avoid Relational errors)
                if not (isinstance(arc_length_simplified, sp.Integral) or arc_length_simplified.has(sp.Integral)):
                    equation = sp.Eq(arc_length_simplified, s_sym)
            except Exception as eq_err:
                try:
                    # Try with original arc_length
                    if not (isinstance(arc_length, sp.Integral) or arc_length.has(sp.Integral)):
                        equation = sp.Eq(arc_length, s_sym)
                except:
                    pass
            
            equation_str = str(equation) if equation else "s = arc_length(t)"
            try:
                equation_latex = sp.latex(equation) if equation else equation_str
            except:
                equation_latex = equation_str
            
            return {
                "msg": "Could not solve for parameter in terms of arc length. The equation may not have an analytical solution. Some curves (like ellipses, certain spirals, or curves with transcendental functions) require numerical methods for arc length reparametrization, which are not currently supported.",
                "equation": equation,
                "equation_str": equation_str,
                "error_type": "NoAnalyticalSolution"
            }
        
        t_in_terms_of_s = solutions[0]
        
        # Validate the solution makes sense
        try:
            # Check if solution is valid (not NaN, not infinite)
            if t_in_terms_of_s == sp.nan or t_in_terms_of_s == sp.oo or t_in_terms_of_s == -sp.oo:
                return {
                    "msg": "Solution found but is invalid (NaN or infinite). The equation may not have a valid analytical solution.",
                    "equation": equation,
                    "error_type": "InvalidSolution"
                }
        except:
            pass

        # Step 3 - Substitute t back into the parametrization
        reparametrized_curve = []
        for coord in parametrization:
            try:
                substituted = coord.subs(parameter, t_in_terms_of_s)
                reparametrized_curve.append(substituted)
            except Exception as e:
                return {
                    "msg": f"Error substituting parameter in parametrization: {str(e)}",
                    "equation": equation,
                }

        # Step 4 - Simplify
        try:
            reparametrized_curve = [sp.simplify(coord) for coord in reparametrized_curve]
        except Exception as e:
            # If simplification fails, return the unsimplified version with a warning
            return {
                "msg": f"Reparametrization computed but simplification failed: {str(e)}",
                "reparametrized_curve": reparametrized_curve,
            }

        return reparametrized_curve
    except Exception as e:
        return {
            "msg": f"Error during reparametrization: {str(e) if str(e) else type(e).__name__}",
            "error_type": type(e).__name__,
        }


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
    # Check if curve is 2D or 3D
    is_2d = len(parametrization) == 2

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
    if is_2d:
        # For 2D curves, use the 2D curvature formula
        # kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
        x_t, y_t = X_t[0], X_t[1]
        x_tt, y_tt = X_tt[0], X_tt[1]
        kappa = sp.Abs(x_t * y_tt - y_t * x_tt) / (X_t_norm**3)
        kappa = kappa.replace(sp.Abs, lambda x: x.args[0] if len(x.args) > 0 else x)
    else:
        # For 3D curves, use cross product
        kappa = (X_t.cross(X_tt)).norm() / (X_t_norm**3)
    kappa = sp.simplify(kappa)

    # Step 5 - Normal vector N
    T_t = Matrix([sp.diff(comp, parameter) for comp in T])
    T_t_norm = sp.sqrt(T_t.dot(T_t))
    if T_t_norm == 0:
        raise ValueError("Normal vector undefined (T' is zero).")
    N = sp.simplify(T_t / T_t_norm)

    # Step 6 - Binormal B
    if is_2d:
        # For 2D curves, binormal is [0, 0, 1] (pointing out of the plane)
        B = Matrix([0, 0, 1])
    else:
        # For 3D curves, compute cross product
        B = sp.simplify(T.cross(N))

    # Step 7 - Torsion tau
    if is_2d:
        # For 2D curves, torsion is always 0 (curve lies in a plane)
        tau = sp.Integer(0)
    else:
        # For 3D curves, compute torsion
        X_ttt = Matrix([sp.diff(coord, parameter) for coord in X_tt])
        cross_product = X_t.cross(X_tt)
        cross_norm_sq = cross_product.norm() ** 2
        if cross_norm_sq == 0:
            # If X_t and X_tt are parallel, torsion is undefined
            tau = sp.Integer(0)
        else:
            tau = cross_product.dot(X_ttt) / cross_norm_sq
        tau = sp.simplify(tau)

    # Pad 2D vectors to 3D for consistency
    if is_2d:
        T = Matrix([T[0], T[1], 0])
        N = Matrix([N[0], N[1], 0])

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
    # Check if the first fundamental form is degenerate (determinant is zero)
    det_I = first_fundamental_form_matrix.det()
    det_I_simplified = sp.simplify(det_I)
    
    # Check if determinant is symbolically zero or can be determined to be zero
    try:
        # Try to check if it's exactly zero
        if det_I_simplified == 0 or sp.simplify(det_I_simplified - 0) == 0:
            raise ValueError("Shape operator undefined for degenerate parametrization: first fundamental form is singular (determinant is zero).")
        # Also check if it simplifies to zero
        if sp.simplify(det_I_simplified).is_zero:
            raise ValueError("Shape operator undefined for degenerate parametrization: first fundamental form is singular (determinant is zero).")
    except (AttributeError, TypeError):
        # If we can't determine symbolically, try numerical check at a sample point
        try:
            # Sample a point to check numerically (u=0, v=0 if possible)
            sample_det = det_I_simplified
            for param in parameters:
                sample_det = sample_det.subs(param, 0)
            if abs(float(sample_det.evalf())) < 1e-10:
                raise ValueError("Shape operator undefined for degenerate parametrization: first fundamental form appears to be singular.")
        except (ValueError, TypeError):
            # If we can't check, proceed but might fail on inversion
            pass
    
    try:
        I_inv = first_fundamental_form_matrix.inv()
        I_inv = sp.simplify(I_inv)
    except Exception as inv_error:
        raise ValueError(f"Shape operator undefined for degenerate parametrization: cannot invert first fundamental form. {str(inv_error)}")

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


def compute_numeric_frenet_serret(t, R):
    import numpy as np

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


def get_curve_display_params(curve, params, t0_str, t1_str):
    exprs = params.get("exprs", {})
    var_provided = params.get("var", "").strip()
    x_provided = exprs.get("x", "").strip() if exprs.get("x") else ""
    y_provided = exprs.get("y", "").strip() if exprs.get("y") else ""
    z_provided = exprs.get("z", "").strip() if exprs.get("z") else ""
    
    default_exprs = get_default_curve_expressions(curve, params)
    
    has_custom_input = (
        (var_provided and var_provided != "t") or
        (x_provided and x_provided != "0") or
        (y_provided and y_provided != "0") or
        (z_provided and z_provided != "0")
    )
    
    if not has_custom_input:
        variable_param = "t"
        x_expr = default_exprs["x"]
        y_expr = default_exprs["y"]
        z_expr = default_exprs["z"]
    else:
        variable_param = var_provided if var_provided else "t"
        x_expr = x_provided if x_provided and x_provided != "0" else default_exprs["x"]
        y_expr = y_provided if y_provided and y_provided != "0" else default_exprs["y"]
        z_expr = z_provided if z_provided and z_provided != "0" else default_exprs["z"]
    
    x_expr = substitute_curve_params_in_expr(x_expr, params, curve)
    y_expr = substitute_curve_params_in_expr(y_expr, params, curve)
    z_expr = substitute_curve_params_in_expr(z_expr, params, curve)
    
    variable_param_str = f'["{variable_param}"]'
    parametrization_str = f'["{x_expr}", "{y_expr}", "{z_expr}"]'
    trange_str = f'["{t0_str}", "{t1_str}"]'
    
    return {
        "variable_param": variable_param_str,
        "parametrization": parametrization_str,
        "trange": trange_str,
    }


def get_computation_steps(computation_type, parametrization=None, parameters=None):
    steps = []
    
    if computation_type == "arc_length":
        steps = [
            "Step 1: Convert the parametrization X(t) = [x(t), y(t), z(t)] into a matrix representation",
            "Step 2: Compute the first derivative X'(t) = [x'(t), y'(t), z'(t)] by differentiating each component with respect to t",
            "Step 3: Calculate the magnitude (speed) of the derivative: ||X'(t)|| = √(x'(t)² + y'(t)² + z'(t)²)",
            "Step 4: Simplify the magnitude expression by removing absolute values and combining like terms",
            "Step 5: Set up the arc length integral: s(t) = ∫||X'(t)|| dt from t₀ to t₁",
            "Step 6: Evaluate the definite integral to obtain the arc length function s(t)",
            "Step 7: Simplify the final arc length expression"
        ]
    
    elif computation_type == "reparam_arc_length":
        steps = [
            "Step 1: Compute the arc length function s(t) = ∫||X'(t)|| dt from t₀ to t₁",
            "Step 2: Simplify the arc length expression s(t) to prepare for solving",
            "Step 3: Set up the equation s(t) = s, where s is the arc length parameter",
            "Step 4: Attempt to solve the equation for t in terms of s using multiple methods:",
            "        - Standard solve method",
            "        - Manual solve with extended options",
            "        - Indefinite integral approach (if speed is constant)",
            "Step 5: Filter solutions to keep only real-valued solutions (remove complex solutions)",
            "Step 6: Select the appropriate solution t(s) from the valid solutions",
            "Step 7: Substitute t(s) back into the original parametrization: X(s) = X(t(s))",
            "Step 8: Simplify each component of the reparametrized curve to get the final result"
        ]
    
    elif computation_type == "frenet":
        steps = [
            "Step 1: Determine if the curve is 2D (two components) or 3D (three components)",
            "Step 2: Compute the first derivative (velocity vector) X'(t) = [x'(t), y'(t), z'(t)]",
            "Step 3: Simplify the velocity vector X'(t)",
            "Step 4: Calculate the magnitude ||X'(t)|| = √(X'(t) · X'(t))",
            "Step 5: Compute the unit tangent vector T = X'(t) / ||X'(t)|| and simplify",
            "Step 6: Compute the second derivative (acceleration) X''(t) = [x''(t), y''(t), z''(t)]",
            "Step 7: Calculate curvature κ:",
            "        - For 2D curves: κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)",
            "        - For 3D curves: κ = ||X'(t) × X''(t)|| / ||X'(t)||³",
            "Step 8: Compute the derivative of the tangent vector T'(t) = dT/dt",
            "Step 9: Calculate the magnitude ||T'(t)|| of the tangent derivative",
            "Step 10: Compute the unit normal vector N = T'(t) / ||T'(t)|| and simplify",
            "Step 11: Calculate the binormal vector B:",
            "         - For 2D curves: B = [0, 0, 1] (pointing out of the plane)",
            "         - For 3D curves: B = T × N (cross product)",
            "Step 12: Compute torsion τ:",
            "         - For 2D curves: τ = 0 (curve lies in a plane)",
            "         - For 3D curves: τ = (X'(t) × X''(t)) · X'''(t) / ||X'(t) × X''(t)||²",
            "Step 13: Pad 2D vectors to 3D format for consistency (add z=0 component)",
            "Step 14: Simplify all components (T, N, B, κ, τ) to get the final Frenet-Serret apparatus"
        ]
    
    elif computation_type == "first_form":
        steps = [
            "Step 1: Extract the parametrization components X(u,v) = [x(u,v), y(u,v), z(u,v)]",
            "Step 2: Compute the partial derivative with respect to u: X_u = ∂X/∂u = [∂x/∂u, ∂y/∂u, ∂z/∂u]",
            "Step 3: Compute the partial derivative with respect to v: X_v = ∂X/∂v = [∂x/∂v, ∂y/∂v, ∂z/∂v]",
            "Step 4: Calculate the first fundamental form coefficient E = X_u · X_u (dot product of X_u with itself)",
            "Step 5: Calculate the first fundamental form coefficient F = X_u · X_v (dot product of X_u and X_v)",
            "Step 6: Calculate the first fundamental form coefficient G = X_v · X_v (dot product of X_v with itself)",
            "Step 7: Simplify each coefficient E, F, and G separately",
            "Step 8: Construct the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 9: Verify the matrix is symmetric (F appears in both off-diagonal positions)"
        ]
    
    elif computation_type == "second_form":
        steps = [
            "Step 1: Compute the first partial derivatives X_u = ∂X/∂u and X_v = ∂X/∂v",
            "Step 2: Calculate the cross product X_u × X_v to get the normal vector direction",
            "Step 3: Compute the magnitude ||X_u × X_v|| of the cross product",
            "Step 4: Check for degenerate parametrization (if ||X_u × X_v|| = 0, the surface is degenerate)",
            "Step 5: Calculate the unit normal vector n = (X_u × X_v) / ||X_u × X_v||",
            "Step 6: Compute the second partial derivatives:",
            "        - X_uu = ∂²X/∂u² = ∂(X_u)/∂u",
            "        - X_uv = ∂²X/∂u∂v = ∂(X_u)/∂v",
            "        - X_vv = ∂²X/∂v² = ∂(X_v)/∂v",
            "Step 7: Calculate the second fundamental form coefficient L = n · X_uu (normal component of X_uu)",
            "Step 8: Calculate the second fundamental form coefficient M = n · X_uv (normal component of X_uv)",
            "Step 9: Calculate the second fundamental form coefficient N = n · X_vv (normal component of X_vv)",
            "Step 10: Simplify each coefficient L, M, and N separately",
            "Step 11: Construct the second fundamental form matrix II = [[L, M], [M, N]]",
            "Step 12: Verify the matrix is symmetric (M appears in both off-diagonal positions)"
        ]
    
    elif computation_type == "gaussian_curvature":
        steps = [
            "Step 1: Compute the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 2: Compute the second fundamental form matrix II = [[L, M], [M, N]]",
            "Step 3: Calculate the determinant of the first fundamental form: det(I) = EG - F²",
            "Step 4: Check for degenerate parametrization (if det(I) = 0, the surface is degenerate)",
            "Step 5: Compute the inverse of the first fundamental form: I⁻¹ = (1/det(I)) · [[G, -F], [-F, E]]",
            "Step 6: Calculate the shape operator matrix S = I⁻¹ · II (matrix multiplication)",
            "Step 7: Compute the Gaussian curvature K = det(S) = det(I⁻¹ · II)",
            "Step 8: Simplify using the formula K = (LN - M²) / (EG - F²)",
            "Step 9: Simplify the final expression for K"
        ]
    
    elif computation_type == "mean_curvature":
        steps = [
            "Step 1: Compute the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 2: Compute the second fundamental form matrix II = [[L, M], [M, N]]",
            "Step 3: Calculate the determinant of the first fundamental form: det(I) = EG - F²",
            "Step 4: Check for degenerate parametrization (if det(I) = 0, the surface is degenerate)",
            "Step 5: Compute the inverse of the first fundamental form: I⁻¹ = (1/det(I)) · [[G, -F], [-F, E]]",
            "Step 6: Calculate the shape operator matrix S = I⁻¹ · II (matrix multiplication)",
            "Step 7: Compute the trace of the shape operator: trace(S) = S₁₁ + S₂₂",
            "Step 8: Calculate the mean curvature H = (1/2) · trace(S)",
            "Step 9: Simplify using the formula H = (EN - 2FM + GL) / (2(EG - F²))",
            "Step 10: Simplify the final expression for H"
        ]
    
    elif computation_type == "principal_curvatures":
        steps = [
            "Step 1: Compute the Gaussian curvature K using the shape operator determinant",
            "Step 2: Compute the mean curvature H using half the trace of the shape operator",
            "Step 3: Calculate the discriminant Δ = √(H² - K)",
            "Step 4: Verify that H² ≥ K (ensures real principal curvatures)",
            "Step 5: Compute the first principal curvature k₁ = H + Δ",
            "Step 6: Compute the second principal curvature k₂ = H - Δ",
            "Step 7: Simplify the expression for k₁",
            "Step 8: Simplify the expression for k₂",
            "Step 9: Verify that k₁ · k₂ = K (product equals Gaussian curvature)",
            "Step 10: Verify that (k₁ + k₂) / 2 = H (average equals mean curvature)"
        ]
    
    elif computation_type == "christoffel":
        steps = [
            "Step 1: Compute the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 2: Extract the metric coefficients E, F, and G from the first fundamental form",
            "Step 3: Compute partial derivatives of the metric coefficients:",
            "        - E_u = ∂E/∂u, E_v = ∂E/∂v",
            "        - F_u = ∂F/∂u, F_v = ∂F/∂v",
            "        - G_u = ∂G/∂u, G_v = ∂G/∂v",
            "Step 4: Calculate the determinant det(I) = EG - F²",
            "Step 5: Compute the inverse of the first fundamental form: I⁻¹ = (1/det(I)) · [[G, -F], [-F, E]]",
            "Step 6: Construct multiplier vectors for each Christoffel symbol:",
            "        - For Γᵤᵤᵘ and Γᵤᵤᵛ: [½E_u, F_u - ½E_v]ᵀ",
            "        - For Γᵤᵥᵘ and Γᵤᵥᵛ: [½E_v, ½G_u]ᵀ",
            "        - For Γᵥᵥᵘ and Γᵥᵥᵛ: [F_v - ½G_u, ½G_v]ᵀ",
            "Step 7: Compute each Christoffel symbol using matrix multiplication:",
            "        - Γᵤᵤᵘ = I⁻¹[0,:] · multiplier_uu",
            "        - Γᵤᵤᵛ = I⁻¹[1,:] · multiplier_uu",
            "        - Γᵤᵥᵘ = I⁻¹[0,:] · multiplier_uv",
            "        - Γᵤᵥᵛ = I⁻¹[1,:] · multiplier_uv",
            "        - Γᵥᵥᵘ = I⁻¹[0,:] · multiplier_vv",
            "        - Γᵥᵥᵛ = I⁻¹[1,:] · multiplier_vv",
            "Step 8: Simplify all six Christoffel symbols (Γᵤᵤᵘ, Γᵤᵤᵛ, Γᵤᵥᵘ, Γᵤᵥᵛ, Γᵥᵥᵘ, Γᵥᵥᵛ)"
        ]
    
    elif computation_type == "gauss_equations":
        steps = [
            "Step 1: Compute the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 2: Compute the second fundamental form matrix II = [[L, M], [M, N]]",
            "Step 3: Compute all six Christoffel symbols (Γᵤᵤᵘ, Γᵤᵤᵛ, Γᵤᵥᵘ, Γᵤᵥᵛ, Γᵥᵥᵘ, Γᵥᵥᵛ)",
            "Step 4: Calculate partial derivatives of the Christoffel symbols with respect to u and v",
            "Step 5: Compute the Riemann curvature tensor components Rᵢⱼₖˡ using:",
            "        - Rᵤᵥᵤᵘ = ∂Γᵤᵥᵘ/∂v - ∂Γᵤᵤᵘ/∂u + Γᵤᵥᵘ·Γᵤᵤᵘ + Γᵤᵥᵛ·Γᵥᵤᵘ - Γᵤᵤᵘ·Γᵤᵥᵘ - Γᵤᵤᵛ·Γᵥᵥᵘ",
            "        - Similar expressions for other components",
            "Step 6: Apply the Gauss equations which relate the Riemann tensor to the second fundamental form:",
            "        - Rᵤᵥᵤᵘ = L·N - M² (first Gauss equation)",
            "        - Rᵤᵥᵤᵛ = ... (second Gauss equation)",
            "        - Rᵤᵥᵥᵘ = ... (third Gauss equation)",
            "        - Rᵤᵥᵥᵛ = ... (fourth Gauss equation)",
            "Step 7: Simplify each of the four Gauss equations",
            "Step 8: Verify the equations are satisfied (left-hand side equals right-hand side)"
        ]
    
    elif computation_type == "codazzi_equations":
        steps = [
            "Step 1: Compute the first fundamental form matrix I = [[E, F], [F, G]]",
            "Step 2: Compute the second fundamental form matrix II = [[L, M], [M, N]]",
            "Step 3: Extract the second fundamental form coefficients L, M, and N",
            "Step 4: Compute all six Christoffel symbols (Γᵤᵤᵘ, Γᵤᵤᵛ, Γᵤᵥᵘ, Γᵤᵥᵛ, Γᵥᵥᵘ, Γᵥᵥᵛ)",
            "Step 5: Calculate partial derivatives of the second fundamental form coefficients:",
            "        - L_v = ∂L/∂v, M_u = ∂M/∂u, M_v = ∂M/∂v, N_u = ∂N/∂u",
            "Step 6: Construct the right-hand side of the first Codazzi equation:",
            "        RHS₁ = L·Γᵤᵥᵘ + M·Γᵤᵥᵛ - M·Γᵤᵤᵘ - N·Γᵤᵤᵛ",
            "Step 7: Construct the right-hand side of the second Codazzi equation:",
            "        RHS₂ = L·Γᵥᵥᵘ + M·Γᵥᵥᵛ - M·Γᵤᵥᵘ - N·Γᵤᵥᵛ",
            "Step 8: Calculate the left-hand side of the first Codazzi equation:",
            "        LHS₁ = L_v - M_u",
            "Step 9: Calculate the left-hand side of the second Codazzi equation:",
            "        LHS₂ = M_v - N_u",
            "Step 10: Simplify both left-hand sides and right-hand sides",
            "Step 11: Verify the Codazzi equations: LHS₁ = RHS₁ and LHS₂ = RHS₂"
        ]
    
    return steps


def get_surface_display_params(surface, params, var_u, var_v, exprs_num, u0_str, u1_str, v0_str, v1_str):
    """
    Get display parameters for surfaces: variable_param, parametrization, urange, and vrange.
    Returns defaults from param_script.js ONLY if no user inputs provided.
    
    Args:
        surface: Surface type name (e.g., "torus", "sphere")
        params: Dictionary of parameters from request
        var_u: Variable name for u parameter (used for mesh generation)
        var_v: Variable name for v parameter (used for mesh generation)
        exprs_num: Dictionary with substituted expressions (from get_surface_expressions)
        u0_str: String representation of u0 range start
        u1_str: String representation of u1 range end
        v0_str: String representation of v0 range start
        v1_str: String representation of v1 range end
    
    Returns:
        dict with keys:
            - variable_param: String list format like '["u", "v"]'
            - parametrization: String list format like '["cos(u)", "sin(u)", "v"]'
            - urange: String list format like '["0", "2*pi"]'
            - vrange: String list format like '["0", "2*pi"]'
    """
    default_surface_exprs = get_default_surface_expressions(surface, params)
    
    u_provided = params.get("u", "").strip()
    v_provided = params.get("v", "").strip()
    x_provided = params.get("x", "").strip()
    y_provided = params.get("y", "").strip()
    z_provided = params.get("z", "").strip()
    
    has_custom_input = (
        (u_provided and u_provided != "u") or
        (v_provided and v_provided != "v") or
        (x_provided and x_provided != "0") or
        (y_provided and y_provided != "0") or
        (z_provided and z_provided != "0")
    )
    
    if not has_custom_input:
        return_var_u = "u"
        return_var_v = "v"
        
        x_expr = substitute_params_in_expr(default_surface_exprs["x"], params, surface)
        y_expr = substitute_params_in_expr(default_surface_exprs["y"], params, surface)
        z_expr = substitute_params_in_expr(default_surface_exprs["z"], params, surface)
    else:
        return_var_u = var_u
        return_var_v = var_v
       
        x_expr = exprs_num.get("x", substitute_params_in_expr(default_surface_exprs["x"], params, surface))
        y_expr = exprs_num.get("y", substitute_params_in_expr(default_surface_exprs["y"], params, surface))
        z_expr = exprs_num.get("z", substitute_params_in_expr(default_surface_exprs["z"], params, surface))
    
    parameters_str = f'["{return_var_u}", "{return_var_v}"]'
    parametrization_str = f'["{x_expr}", "{y_expr}", "{z_expr}"]'
    urange_str = f'["{u0_str}", "{u1_str}"]'
    vrange_str = f'["{v0_str}", "{v1_str}"]'
    
    return {
        "variable_param": parameters_str,
        "parametrization": parametrization_str,
        "urange": urange_str,
        "vrange": vrange_str,
    }
