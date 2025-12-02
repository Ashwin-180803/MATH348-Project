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

    # CURVE
    if mode == "curve":
        try:
            curve = data["curve"]
            params = data.get("params", {})
            t0 = float(data.get("t0", 0.0))
            t1 = float(data.get("t1", 2 * 3.141592653589793))
            n = int(data.get("n", 400))
            symbolic_quantity = data.get("symbolic_quantity", "")

            # Numeric: positions, curvature, torsion, arc-length
            t, R = functions.numeric_curve_positions(curve, params, t0, t1, n)
            curvature, torsion = functions.curvature_torsion_from_R(R, t)
            arclen = functions.arc_length_numeric(R, t)

            symbolic = {}

            if symbolic_quantity:
                t_sym = sp.symbols("t", real=True)
                r_matrix = functions.symbolic_formula_for(curve, params)
                if isinstance(r_matrix, sp.Matrix):
                    param_list = list(r_matrix)

                    # Parse input
                    for i in range(len(param_list)):
                        param_list[i] = parser.parse_input(str(param_list[i]))

                else:
                    # fallback just in case
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

            return JsonResponse(
                {
                    "ok": True,
                    "mode": "curve",
                    "t": t.tolist(),
                    "x": R[:, 0].tolist(),
                    "y": R[:, 1].tolist(),
                    "z": R[:, 2].tolist(),
                    "curvature": curvature,
                    "torsion": torsion,
                    "arc_length": arclen,
                    "symbolic": symbolic,
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
            u0 = float(params.get("u0", 0.0))
            u1 = float(params.get("u1", 2 * 3.141592653589793))
            v0 = float(params.get("v0", 0.0))
            v1 = float(params.get("v1", 2 * 3.141592653589793))

            exprs_num = functions.get_surface_expressions(surface, params)
            
            U, V, X, Y, Z = functions.mesh_from_parametric_surfaces(
                exprs_num, (u0, u1), (v0, v1), nu, nv
            )
            
            symbolic = {}

            if compute_symbolic_flag and symbolic_quantity:
                u, v = sp.symbols("u v", real=True)

                param_list = [
                    sp.sympify(exprs_num["x"]),
                    sp.sympify(exprs_num["y"]),
                    sp.sympify(exprs_num["z"]),
                ]
                
                parameters = (u, v)
                # Parse input
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

            return JsonResponse(
                {
                    "ok": True,
                    "mode": "surface",
                    "X": X.tolist(),
                    "Y": Y.tolist(),
                    "Z": Z.tolist(),
                    "symbolic": symbolic,
                }
            )

        except Exception as e:
            return JsonResponse({"ok": False, "error": str(e)}, status=400)

    return JsonResponse({"ok": False, "error": "Unknown mode"}, status=400)
