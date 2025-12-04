import json
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from calcapp.utils import functions, parser
import sympy as sp


def index(request):
    return render(request, "index.html")


@csrf_exempt
def compute(request):
    if request.method != "POST":
        return JsonResponse({"ok": False, "error": "POST only"}, status=400)

    try:
        data = json.loads(request.body.decode("utf-8"))
    except Exception:
        return JsonResponse({"ok": False, "error": "Invalid JSON"}, status=400)

    mode = data.get("mode", "curve")
    print(data)
    # CURVE
    if mode == "curve":
        try:
            curve = data["curve"]
            params = data.get("params", {})
            t0_str = data.get("t0", "0")
            t1_str = data.get("t1", "2*pi")
            t0 = functions.parse_range_value(t0_str, 0.0)
            t1 = functions.parse_range_value(t1_str, 2 * 3.141592653589793)
            n = int(data.get("n", 400))
            symbolic_quantity = data.get("symbolic_quantity", "")

            
            t, R = functions.numeric_curve_positions(curve, params, t0, t1, n)

            
            frenet_data = functions.compute_numeric_frenet_serret(t, R)

            symbolic = {}

            if symbolic_quantity:
                t_sym = sp.symbols("t", real=True)
                r_matrix = functions.symbolic_formula_for(curve, params)
                if isinstance(r_matrix, sp.Matrix):
                    param_list = list(r_matrix)

                    
                    for i in range(len(param_list)):
                        param_list[i] = parser.parse_input(str(param_list[i]))

                else:
                    
                    param_list = list(sp.Matrix(r_matrix))

                bounds = {t_sym: (t0, t1)}

                if symbolic_quantity == "arc_length":
                    try:
                        s_expr = functions.compute_arc_length(param_list, t_sym, bounds)
                        symbolic["s(t)"] = sp.latex(s_expr)
                    except Exception as e:
                        symbolic["arc_length_error"] = str(e)

                elif symbolic_quantity == "reparam_arc_length":
                    try:
                        rep = functions.compute_arc_length_reparametrization(
                            param_list, t_sym, bounds
                        )

                        s = sp.symbols("s", real=True)
                        for i, coord in enumerate(rep):
                            symbolic[f"x_{i + 1}(s)"] = sp.latex(coord)
                    except Exception as e:
                        symbolic["reparam_error"] = str(e)

                elif symbolic_quantity == "frenet":
                    try:
                        fr = functions.compute_frenet_serret_apparatus(
                            param_list, t_sym
                        )
                        for k, v in fr.items():
                            symbolic[k] = sp.latex(v)
                    except Exception as e:
                        symbolic["frenet_error"] = str(e)

            
            exprs = params.get("exprs", {})
            var_provided = params.get("var", "").strip()
            x_provided = exprs.get("x", "").strip() if exprs.get("x") else ""
            y_provided = exprs.get("y", "").strip() if exprs.get("y") else ""
            z_provided = exprs.get("z", "").strip() if exprs.get("z") else ""
            
            
            default_exprs = functions.get_default_curve_expressions(curve, params)
            
            
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
            
            
            variable_param_str = f'["{variable_param}"]'
            parametrization_str = f'["{x_expr}", "{y_expr}", "{z_expr}"]'
            trange_str = f'["{t0_str}", "{t1_str}"]'
            
            return JsonResponse(
                {
                    "ok": True,
                    "mode": "curve",
                    "t": t.tolist(),
                    "x": R[:, 0].tolist(),
                    "y": R[:, 1].tolist(),
                    "z": R[:, 2].tolist(),
                    "curvature": frenet_data["curvature"],
                    "torsion": frenet_data["torsion"],
                    "arc_length": frenet_data["arc_length"],
                    "symbolic": symbolic,
                    "variable_param": variable_param_str,
                    "parametrization": parametrization_str,
                    "trange": trange_str,
                }
            )
        except Exception as e:
            return JsonResponse({"ok": False, "error": str(e)}, status=400)

    # SURFACE
    if mode == "surface":
        try:
            surface = data["surface"]
            params = data.get("params", {})
            compute_symbolic_flag = bool(data.get("compute_symbolic", False))
            symbolic_quantity = data.get("symbolic_quantity", "")

            nu = int(params.get("nu", 60))
            nv = int(params.get("nv", 60))
            u0_str = params.get("u0", "0")
            u1_str = params.get("u1", "2*pi")
            v0_str = params.get("v0", "0")
            v1_str = params.get("v1", "2*pi")
            u0 = functions.parse_range_value(u0_str, 0.0)
            u1 = functions.parse_range_value(u1_str, 2 * 3.141592653589793)
            v0 = functions.parse_range_value(v0_str, 0.0)
            v1 = functions.parse_range_value(v1_str, 2 * 3.141592653589793)

            # Get variable names (will be updated later if using defaults)
            var_u = params.get("u", "u")
            var_v = params.get("v", "v")

            exprs_num = functions.get_surface_expressions(surface, params)
            print(f"DEBUG: Surface {surface}, variables: u='{var_u}', v='{var_v}', expressions: {exprs_num}")

            U, V, X, Y, Z = functions.mesh_from_parametric_surfaces(
                exprs_num, (u0, u1), (v0, v1), nu, nv, var_u, var_v
            )
            print(f"DEBUG: Mesh generated, X.shape: {X.shape}, Y.shape: {Y.shape}, Z.shape: {Z.shape}")

            
            if X.size == 0 or Y.size == 0 or Z.size == 0:
                raise ValueError("Generated mesh is empty")

            symbolic = {}

            if compute_symbolic_flag and symbolic_quantity:
                u_sym, v_sym = sp.symbols(f"{var_u} {var_v}", real=True)

                param_list = [
                    sp.sympify(exprs_num["x"]),
                    sp.sympify(exprs_num["y"]),
                    sp.sympify(exprs_num["z"]),
                ]

                parameters = (u_sym, v_sym)
                
                for i in range(len(param_list)):
                    param_list[i] = parser.parse_input(str(param_list[i]))

                if symbolic_quantity == "first_form":
                    try:
                        I = functions.compute_first_fundamental_form(
                            param_list, parameters
                        )

                        E = I[0, 0]
                        F = I[0, 1]
                        G = I[1, 1]
                        symbolic["E"] = sp.latex(E)
                        symbolic["F"] = sp.latex(F)
                        symbolic["G"] = sp.latex(G)
                        symbolic["I"] = sp.latex(I)
                    except Exception as e:
                        symbolic["first_form_error"] = str(e)

                elif symbolic_quantity == "second_form":
                    try:
                        II = functions.compute_second_fundamental_form(
                            param_list, parameters
                        )
                        L = II[0, 0]
                        M = II[0, 1]
                        N = II[1, 1]
                        symbolic["L"] = sp.latex(L)
                        symbolic["M"] = sp.latex(M)
                        symbolic["N"] = sp.latex(N)
                        symbolic["II"] = sp.latex(II)
                    except Exception as e:
                        symbolic["second_form_error"] = str(e)

                elif symbolic_quantity in (
                    "gaussian_curvature",
                    "mean_curvature",
                    "principal_curvatures",
                ):
                    try:
                        K_expr = functions.compute_gaussian_curvature(
                            param_list, parameters
                        )
                        H_expr = functions.compute_mean_curvature(
                            param_list, parameters
                        )
                        if symbolic_quantity == "gaussian_curvature":
                            symbolic["K"] = sp.latex(K_expr)
                        elif symbolic_quantity == "mean_curvature":
                            symbolic["H"] = sp.latex(H_expr)
                        else:
                            disc = sp.sqrt(H_expr**2 - K_expr)
                            k1 = sp.simplify(H_expr + disc)
                            k2 = sp.simplify(H_expr - disc)
                            symbolic["H"] = sp.latex(H_expr)
                            symbolic["K"] = sp.latex(K_expr)
                            symbolic["k1"] = sp.latex(k1)
                            symbolic["k2"] = sp.latex(k2)
                    except Exception as e:
                        symbolic["curvature_error"] = str(e)

                elif symbolic_quantity == "christoffel":
                    try:
                        gamma = functions.compute_christoffel_symbols(
                            param_list, parameters
                        )

                        for name, expr in gamma.items():
                            symbolic[name] = sp.latex(expr)
                    except Exception as e:
                        symbolic["christoffel_error"] = str(e)

                elif symbolic_quantity == "gauss_equations":
                    try:
                        g_eqs = functions.compute_gauss_equations(
                            param_list, parameters
                        )
                        for i, expr in enumerate(g_eqs, start=1):
                            symbolic[f"Gauss_eq_{i}"] = sp.latex(expr)
                    except Exception as e:
                        symbolic["gauss_eq_error"] = str(e)

                elif symbolic_quantity == "codazzi_equations":
                    try:
                        c_eqs = functions.compute_codazzi_equations(
                            param_list, parameters
                        )
                        for i, expr in enumerate(c_eqs, start=1):
                            symbolic[f"Codazzi_eq_{i}"] = sp.latex(expr)
                    except Exception as e:
                        symbolic["codazzi_eq_error"] = str(e)

            
            
            default_surface_exprs = functions.get_default_surface_expressions(surface, params)
            
            
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
                
                x_expr = default_surface_exprs["x"]
                y_expr = default_surface_exprs["y"]
                z_expr = default_surface_exprs["z"]
            else:
                
                return_var_u = var_u
                return_var_v = var_v
                x_expr = exprs_num.get("x", default_surface_exprs["x"])
                y_expr = exprs_num.get("y", default_surface_exprs["y"])
                z_expr = exprs_num.get("z", default_surface_exprs["z"])
            
            
            parameters_str = f'["{return_var_u}", "{return_var_v}"]'
            parametrization_str = f'["{x_expr}", "{y_expr}", "{z_expr}"]'
            urange_str = f'["{u0_str}", "{u1_str}"]'
            vrange_str = f'["{v0_str}", "{v1_str}"]'
            
            return JsonResponse(
                {
                    "ok": True,
                    "mode": "surface",
                    "X": X.tolist(),
                    "Y": Y.tolist(),
                    "Z": Z.tolist(),
                    "symbolic": symbolic,
                    "variable_param": parameters_str,
                    "parametrization": parametrization_str,
                    "urange": urange_str,
                    "vrange": vrange_str,
                }
            )

        except Exception as e:
            return JsonResponse({"ok": False, "error": str(e)}, status=400)

    return JsonResponse({"ok": False, "error": "Unknown mode"}, status=400)
