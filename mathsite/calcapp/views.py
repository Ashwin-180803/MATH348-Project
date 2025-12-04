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
            curve = data.get("curve", "")
            if not curve:
                return JsonResponse(
                    {"ok": False, "error": "Curve type not specified"}, status=400
                )

            params = data.get("params", {})
            if not isinstance(params, dict):
                params = {}

            t0_str = data.get("t0", "0")
            t1_str = data.get("t1", "2*pi")
            try:
                t0 = functions.parse_range_value(t0_str, 0.0)
                t1 = functions.parse_range_value(t1_str, 2 * 3.141592653589793)
            except Exception as e:
                return JsonResponse(
                    {"ok": False, "error": f"Error parsing parameter range: {str(e)}"},
                    status=400,
                )

            try:
                n = int(data.get("n", 400))
                if n <= 0:
                    n = 400
                if n > 10000:
                    n = 10000
            except (ValueError, TypeError):
                n = 400

            symbolic_quantity = data.get("symbolic_quantity", "")

            try:
                t, R = functions.numeric_curve_positions(curve, params, t0, t1, n)
            except Exception as e:
                error_msg = (
                    str(e)
                    if str(e)
                    else f"Error computing curve positions: {type(e).__name__}"
                )
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

            try:
                frenet_data = functions.compute_numeric_frenet_serret(t, R)
            except Exception as e:
                error_msg = (
                    str(e)
                    if str(e)
                    else f"Error computing Frenet-Serret data: {type(e).__name__}"
                )
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

            symbolic = {}

            if symbolic_quantity:
                try:
                    t_sym = sp.symbols("t", real=True)
                except Exception as e:
                    symbolic["symbolic_error"] = f"Error creating symbol: {str(e)}"
                    symbolic_quantity = ""

                if symbolic_quantity:
                    try:
                        r_matrix = functions.symbolic_formula_for(curve, params)

                    except Exception as e:
                        symbolic["symbolic_error"] = (
                            f"Error getting symbolic formula: {str(e)}"
                        )
                        symbolic_quantity = ""

                original_symbols = set()
                if isinstance(r_matrix, sp.Matrix):
                    for expr in r_matrix:
                        if hasattr(expr, "free_symbols"):
                            original_symbols.update(expr.free_symbols)
                else:
                    for expr in list(sp.Matrix(r_matrix)):
                        if hasattr(expr, "free_symbols"):
                            original_symbols.update(expr.free_symbols)

                t_symbol_in_expr = None
                for sym in original_symbols:
                    if str(sym) == "t":
                        t_symbol_in_expr = sym
                        break
                if t_symbol_in_expr:
                    original_symbols.discard(t_symbol_in_expr)

                subs_dict_original = {}
                for key, value in params.items():
                    if key != "exprs" and key != "var":
                        try:
                            param_value = float(value)

                            matching_symbol = None
                            for sym in original_symbols:
                                if str(sym) == key:
                                    matching_symbol = sym
                                    break

                            if matching_symbol:
                                subs_dict_original[matching_symbol] = param_value
                        except (ValueError, TypeError):
                            pass

                if isinstance(r_matrix, sp.Matrix):
                    param_list = list(r_matrix)
                else:
                    param_list = list(sp.Matrix(r_matrix))

                if subs_dict_original:
                    for i in range(len(param_list)):
                        param_list[i] = param_list[i].subs(subs_dict_original)

                        for key, value in params.items():
                            if key != "exprs" and key != "var":
                                try:
                                    param_value = float(value)

                                    for sym in param_list[i].free_symbols:
                                        if str(sym) == key:
                                            param_list[i] = param_list[i].xreplace(
                                                {sym: param_value}
                                            )
                                except (ValueError, TypeError, AttributeError):
                                    pass

                        param_list[i] = sp.simplify(param_list[i])

                try:
                    for i in range(len(param_list)):
                        try:
                            param_list[i] = parser.parse_input(
                                str(param_list[i]), [t_sym]
                            )
                        except Exception as parse_err:
                            symbolic["parse_error"] = (
                                f"Error parsing expression {i + 1}: {str(parse_err)}"
                            )
                            symbolic_quantity = ""  # Skip symbolic computation
                            break
                except Exception as e:
                    symbolic["parse_error"] = f"Error during parsing: {str(e)}"
                    symbolic_quantity = ""  # Skip symbolic computation

                all_symbols_after = set()
                for expr in param_list:
                    if hasattr(expr, "free_symbols"):
                        all_symbols_after.update(expr.free_symbols)
                all_symbols_after.discard(t_sym)

                subs_dict_after = {}
                for key, value in params.items():
                    if key != "exprs" and key != "var":
                        try:
                            param_value = float(value)

                            for sym in all_symbols_after:
                                if str(sym) == key:
                                    subs_dict_after[sym] = param_value
                        except (ValueError, TypeError):
                            pass

                if subs_dict_after:
                    for i in range(len(param_list)):
                        param_list[i] = param_list[i].subs(subs_dict_after)
                        param_list[i] = sp.simplify(param_list[i])

                bounds = [t0, t1]

                if symbolic_quantity == "arc_length":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "arc_length", param_list, [t_sym]
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        s_expr = functions.compute_arc_length(param_list, t_sym, bounds)
                        print(s_expr)
                        if isinstance(s_expr, dict):
                            symbolic["arc_length_error"] = s_expr.get(
                                "msg", "No elementary antiderivative found."
                            )
                            try:
                                symbolic["arc_length_integral"] = sp.latex(
                                    s_expr.get("integral", "")
                                )
                            except:
                                symbolic["arc_length_integral"] = str(
                                    s_expr.get("integral", "")
                                )
                        else:
                            try:
                                symbolic["s(t)"] = sp.latex(s_expr)
                            except Exception as latex_err:
                                symbolic["s(t)"] = str(s_expr)
                                symbolic["arc_length_warning"] = (
                                    f"LaTeX conversion failed: {str(latex_err)}"
                                )
                    except Exception as e:
                        symbolic["arc_length_error"] = (
                            f"Error computing arc length: {str(e)}"
                        )

                elif symbolic_quantity == "reparam_arc_length":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "reparam_arc_length", param_list, [t_sym]
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass
                        rep = functions.compute_arc_length_reparametrization(
                            param_list, t_sym, bounds
                        )
                        if isinstance(rep, dict):
                            error_msg = rep.get("msg", "Reparametrization failed.")
                            if not error_msg or error_msg.strip() == "":
                                error_msg = f"Reparametrization failed. Error type: {rep.get('error_type', 'Unknown')}"

                            if "analytical solution" in error_msg.lower():
                                error_msg = "Arc length reparametrization failed: This curve does not have a closed-form analytical solution. Some curves (like ellipses, certain spirals, or curves with transcendental functions) require numerical methods for arc length reparametrization, which are not currently supported."
                            elif (
                                "No elementary antiderivative" in error_msg
                                or "antiderivative" in error_msg.lower()
                            ):
                                error_msg = "Arc length reparametrization failed: The arc length integral for this curve does not have an elementary antiderivative. This means the arc length cannot be expressed in closed form using standard functions. Curves like ellipses, certain polynomials of degree 4 or higher, and curves with complex expressions often have this limitation."

                            symbolic["reparam_error"] = error_msg

                            if "integral" in rep:
                                try:
                                    symbolic["reparam_integral"] = sp.latex(
                                        rep["integral"]
                                    )
                                    symbolic["reparam_integral_note"] = (
                                        "The arc length is given by this integral, which cannot be evaluated in closed form:"
                                    )
                                except:
                                    symbolic["reparam_integral"] = str(rep["integral"])

                            if "equation" in rep:
                                try:
                                    symbolic["reparam_equation"] = sp.latex(
                                        rep["equation"]
                                    )
                                except:
                                    if "equation_str" in rep:
                                        symbolic["reparam_equation"] = rep[
                                            "equation_str"
                                        ]
                                    else:
                                        symbolic["reparam_equation"] = str(
                                            rep["equation"]
                                        )

                            if "reparametrized_curve" in rep:
                                s = sp.symbols("s", real=True)
                                for i, coord in enumerate(rep["reparametrized_curve"]):
                                    try:
                                        symbolic[f"x_{i + 1}(s)"] = sp.latex(coord)
                                    except:
                                        symbolic[f"x_{i + 1}(s)"] = str(coord)
                        else:
                            s = sp.symbols("s", real=True)
                            for i, coord in enumerate(rep):
                                try:
                                    symbolic[f"x_{i + 1}(s)"] = sp.latex(coord)
                                except Exception as latex_err:
                                    symbolic[f"x_{i + 1}(s)"] = str(coord)
                                    symbolic[f"reparam_warning_{i + 1}"] = (
                                        f"LaTeX conversion failed: {str(latex_err)}"
                                    )
                    except Exception as e:
                        error_msg = (
                            str(e) if str(e) else f"{type(e).__name__}: {repr(e)}"
                        )
                        symbolic["reparam_error"] = (
                            f"Error computing reparametrization: {error_msg}"
                        )

                elif symbolic_quantity == "frenet":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "frenet", param_list, [t_sym]
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass
                        fr = functions.compute_frenet_serret_apparatus(
                            param_list, t_sym
                        )

                        for k, v in fr.items():
                            try:
                                symbolic[k] = sp.latex(v)
                            except Exception as latex_err:
                                symbolic[k] = str(v)
                                symbolic[f"{k}_latex_warning"] = (
                                    f"LaTeX conversion failed: {str(latex_err)}"
                                )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["frenet_error"] = (
                            f"Error computing Frenet-Serret apparatus: {error_msg}"
                        )

            try:
                display_params = functions.get_curve_display_params(
                    curve, params, t0_str, t1_str
                )
            except Exception as e:
                display_params = {
                    "variable_param": '["t"]',
                    "parametrization": '["", "", ""]',
                    "trange": f'["{t0_str}", "{t1_str}"]',
                }

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
                    "variable_param": display_params["variable_param"],
                    "parametrization": display_params["parametrization"],
                    "trange": display_params["trange"],
                }
            )
        except Exception as e:
            error_msg = f"Error serializing response: {str(e)}"
            return JsonResponse({"ok": False, "error": error_msg}, status=400)
        except Exception as e:
            return JsonResponse({"ok": False, "error": str(e)}, status=400)

    # SURFACE
    if mode == "surface":
        try:
            surface = data.get("surface", "")
            if not surface:
                return JsonResponse(
                    {"ok": False, "error": "Surface type not specified"}, status=400
                )

            params = data.get("params", {})
            if not isinstance(params, dict):
                params = {}

            compute_symbolic_flag = bool(data.get("compute_symbolic", False))
            symbolic_quantity = data.get("symbolic_quantity", "")

            try:
                nu = int(params.get("nu", 60))
                if nu <= 0:
                    nu = 60
                if nu > 200:
                    nu = 200
            except (ValueError, TypeError):
                nu = 60

            try:
                nv = int(params.get("nv", 60))
                if nv <= 0:
                    nv = 60
                if nv > 200:
                    nv = 200
            except (ValueError, TypeError):
                nv = 60

            u0_str = params.get("u0", "0")
            u1_str = params.get("u1", "2*pi")
            v0_str = params.get("v0", "0")
            v1_str = params.get("v1", "2*pi")

            try:
                u0 = functions.parse_range_value(u0_str, 0.0)
                u1 = functions.parse_range_value(u1_str, 2 * 3.141592653589793)
                v0 = functions.parse_range_value(v0_str, 0.0)
                v1 = functions.parse_range_value(v1_str, 2 * 3.141592653589793)
            except Exception as e:
                return JsonResponse(
                    {"ok": False, "error": f"Error parsing parameter range: {str(e)}"},
                    status=400,
                )

            var_u = params.get("u", "u")
            var_v = params.get("v", "v")
            if not isinstance(var_u, str) or not isinstance(var_v, str):
                var_u = "u"
                var_v = "v"

            try:
                exprs_num = functions.get_surface_expressions(surface, params)
            except Exception as e:
                error_msg = (
                    str(e)
                    if str(e)
                    else f"Error getting surface expressions: {type(e).__name__}"
                )
                return JsonResponse({"ok": False, "error": error_msg}, status=400)
            print(
                f"DEBUG: Surface {surface}, variables: u='{var_u}', v='{var_v}', expressions: {exprs_num}"
            )

            try:
                U, V, X, Y, Z = functions.mesh_from_parametric_surfaces(
                    exprs_num, (u0, u1), (v0, v1), nu, nv, var_u, var_v
                )
                print(
                    f"DEBUG: Mesh generated, X.shape: {X.shape}, Y.shape: {Y.shape}, Z.shape: {Z.shape}"
                )
            except Exception as mesh_error:
                error_msg = f"Error generating surface mesh: {str(mesh_error)}"
                print(f"DEBUG: {error_msg}")
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

            if X.size == 0 or Y.size == 0 or Z.size == 0:
                error_msg = f"Generated mesh is empty: X.size={X.size}, Y.size={Y.size}, Z.size={Z.size}"
                print(f"DEBUG: {error_msg}")
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

            try:
                X_list = X.tolist()
                Y_list = Y.tolist()
                Z_list = Z.tolist()

                if not X_list or not Y_list or not Z_list:
                    error_msg = "Arrays converted to empty lists"
                    print(f"DEBUG: {error_msg}")
                    return JsonResponse({"ok": False, "error": error_msg}, status=400)
            except Exception as convert_error:
                error_msg = f"Error converting arrays to lists: {str(convert_error)}"
                print(f"DEBUG: {error_msg}")
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

            symbolic = {}

            if compute_symbolic_flag and symbolic_quantity:
                try:
                    u_sym, v_sym = sp.symbols(f"{var_u} {var_v}", real=True)
                except Exception as e:
                    symbolic["symbolic_error"] = f"Error creating symbols: {str(e)}"
                    compute_symbolic_flag = False

                if compute_symbolic_flag:
                    try:
                        param_list = [
                            sp.sympify(exprs_num.get("x", "0")),
                            sp.sympify(exprs_num.get("y", "0")),
                            sp.sympify(exprs_num.get("z", "0")),
                        ]
                    except Exception as e:
                        symbolic["symbolic_error"] = (
                            f"Error converting expressions: {str(e)}"
                        )
                        compute_symbolic_flag = False

                    if compute_symbolic_flag:
                        parameters = (u_sym, v_sym)

                        try:
                            for i in range(len(param_list)):
                                try:
                                    param_list[i] = parser.parse_input(
                                        str(param_list[i]), list(parameters)
                                    )
                                except Exception as parse_err:
                                    symbolic["parse_error"] = (
                                        f"Error parsing expression {i + 1}: {str(parse_err)}"
                                    )
                                    compute_symbolic_flag = False
                                    break
                        except Exception as e:
                            symbolic["parse_error"] = f"Error during parsing: {str(e)}"
                            compute_symbolic_flag = False

                if compute_symbolic_flag and symbolic_quantity == "first_form":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "first_form", param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        I = functions.compute_first_fundamental_form(
                            param_list, parameters
                        )
                        E = I[0][0, 0]
                        F = I[0][0, 1]
                        G = I[0][1, 1]
                        try:
                            symbolic["E"] = sp.latex(E)
                            symbolic["F"] = sp.latex(F)
                            symbolic["G"] = sp.latex(G)
                            symbolic["I"] = sp.latex(I)
                        except Exception as latex_err:
                            symbolic["E"] = str(E)
                            symbolic["F"] = str(F)
                            symbolic["G"] = str(G)
                            symbolic["I"] = str(I)
                            symbolic["first_form_latex_warning"] = (
                                f"LaTeX conversion failed: {str(latex_err)}"
                            )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["first_form_error"] = (
                            f"Error computing first fundamental form: {error_msg}"
                        )

                elif symbolic_quantity == "second_form":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "second_form", param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        II = functions.compute_second_fundamental_form(
                            param_list, parameters
                        )
                        L = II[0][0, 0]
                        M = II[0][0, 1]
                        N = II[0][1, 1]
                        try:
                            symbolic["L"] = sp.latex(L)
                            symbolic["M"] = sp.latex(M)
                            symbolic["N"] = sp.latex(N)
                            symbolic["II"] = sp.latex(II)
                        except Exception as latex_err:
                            symbolic["L"] = str(L)
                            symbolic["M"] = str(M)
                            symbolic["N"] = str(N)
                            symbolic["II"] = str(II)
                            symbolic["second_form_latex_warning"] = (
                                f"LaTeX conversion failed: {str(latex_err)}"
                            )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["second_form_error"] = (
                            f"Error computing second fundamental form: {error_msg}"
                        )

                elif symbolic_quantity in (
                    "gaussian_curvature",
                    "mean_curvature",
                    "principal_curvatures",
                ):
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                symbolic_quantity, param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        K_expr = functions.compute_gaussian_curvature(
                            param_list, parameters
                        )
                        H_expr = functions.compute_mean_curvature(
                            param_list, parameters
                        )
                        try:
                            if symbolic_quantity == "gaussian_curvature":
                                symbolic["K"] = sp.latex(K_expr)
                            elif symbolic_quantity == "mean_curvature":
                                symbolic["H"] = sp.latex(H_expr)
                            else:
                                disc = sp.sqrt(H_expr**2 - K_expr)
                                k1 = sp.simplify(H_expr + disc)
                                k2 = sp.simplify(H_expr - disc)
                                try:
                                    symbolic["H"] = sp.latex(H_expr)
                                    symbolic["K"] = sp.latex(K_expr)
                                    symbolic["k1"] = sp.latex(k1)
                                    symbolic["k2"] = sp.latex(k2)
                                except Exception as latex_err:
                                    symbolic["H"] = str(H_expr)
                                    symbolic["K"] = str(K_expr)
                                    symbolic["k1"] = str(k1)
                                    symbolic["k2"] = str(k2)
                                    symbolic["curvature_latex_warning"] = (
                                        f"LaTeX conversion failed: {str(latex_err)}"
                                    )
                        except Exception as latex_err:
                            symbolic["curvature_latex_error"] = (
                                f"LaTeX conversion failed: {str(latex_err)}"
                            )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"

                        if (
                            "degenerate" in error_msg.lower()
                            or "singular" in error_msg.lower()
                        ):
                            error_msg = "Curvature cannot be computed: the surface parametrization is degenerate (the first fundamental form is singular). This typically occurs when the partial derivatives are linearly dependent or the surface has zero area at some points."
                        else:
                            error_msg = f"Error computing curvature: {error_msg}"
                        symbolic["curvature_error"] = error_msg

                elif symbolic_quantity == "christoffel":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "christoffel", param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        gamma = functions.compute_christoffel_symbols(
                            param_list, parameters
                        )[0]

                        for name, expr in gamma.items():
                            try:
                                symbolic[name] = sp.latex(expr)
                            except Exception as latex_err:
                                symbolic[name] = str(expr)
                                symbolic[f"{name}_latex_warning"] = (
                                    f"LaTeX conversion failed: {str(latex_err)}"
                                )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["christoffel_error"] = (
                            f"Error computing Christoffel symbols: {error_msg}"
                        )

                elif symbolic_quantity == "gauss_equations":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "gauss_equations", param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        g_eqs = functions.compute_gauss_equations(
                            param_list, parameters
                        )
                        for i, expr in enumerate(g_eqs, start=1):
                            try:
                                symbolic[f"Gauss_eq_{i}"] = sp.latex(expr)
                            except Exception as latex_err:
                                symbolic[f"Gauss_eq_{i}"] = str(expr)
                                symbolic[f"Gauss_eq_{i}_latex_warning"] = (
                                    f"LaTeX conversion failed: {str(latex_err)}"
                                )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["gauss_eq_error"] = (
                            f"Error computing Gauss equations: {error_msg}"
                        )

                elif symbolic_quantity == "codazzi_equations":
                    try:
                        # Add computation steps
                        try:
                            steps = functions.get_computation_steps(
                                "codazzi_equations", param_list, parameters
                            )
                            symbolic["computation_steps"] = steps
                        except Exception:
                            pass

                        c_eqs = functions.compute_codazzi_equations(
                            param_list, parameters
                        )
                        for i, expr in enumerate(c_eqs, start=1):
                            try:
                                symbolic[f"Codazzi_eq_{i}"] = sp.latex(expr)
                            except Exception as latex_err:
                                symbolic[f"Codazzi_eq_{i}"] = str(expr)
                                symbolic[f"Codazzi_eq_{i}_latex_warning"] = (
                                    f"LaTeX conversion failed: {str(latex_err)}"
                                )
                    except Exception as e:
                        error_msg = str(e) if str(e) else f"{type(e).__name__}"
                        symbolic["codazzi_eq_error"] = (
                            f"Error computing Codazzi equations: {error_msg}"
                        )

            default_surface_exprs = functions.get_default_surface_expressions(
                surface, params
            )

            u_provided = params.get("u", "").strip()
            v_provided = params.get("v", "").strip()
            x_provided = params.get("x", "").strip()
            y_provided = params.get("y", "").strip()
            z_provided = params.get("z", "").strip()

            has_custom_input = (
                (u_provided and u_provided != "u")
                or (v_provided and v_provided != "v")
                or (x_provided and x_provided != "0")
                or (y_provided and y_provided != "0")
                or (z_provided and z_provided != "0")
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

            try:
                display_params = functions.get_surface_display_params(
                    surface,
                    params,
                    var_u,
                    var_v,
                    exprs_num,
                    u0_str,
                    u1_str,
                    v0_str,
                    v1_str,
                )
            except Exception as e:
                display_params = {
                    "variable_param": f'["{var_u}", "{var_v}"]',
                    "parametrization": '["", "", ""]',
                    "urange": f'["{u0_str}", "{u1_str}"]',
                    "vrange": f'["{v0_str}", "{v1_str}"]',
                }

            try:
                return JsonResponse(
                    {
                        "ok": True,
                        "mode": "surface",
                        "X": X_list,
                        "Y": Y_list,
                        "Z": Z_list,
                        "symbolic": symbolic,
                        "variable_param": display_params["variable_param"],
                        "parametrization": display_params["parametrization"],
                        "urange": display_params["urange"],
                        "vrange": display_params["vrange"],
                    }
                )
            except Exception as e:
                error_msg = f"Error serializing response: {str(e)}"
                return JsonResponse({"ok": False, "error": error_msg}, status=400)

        except Exception as e:
            return JsonResponse({"ok": False, "error": str(e)}, status=400)

    return JsonResponse({"ok": False, "error": "Unknown mode"}, status=400)
