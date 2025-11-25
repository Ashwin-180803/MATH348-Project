from flask import Flask, request, jsonify, send_from_directory
import numpy as np
import sympy as sp
import math
import os

app = Flask(__name__, static_folder='static', static_url_path='/static')


def numeric_curve_positions(curve_name, params, t0, t1, n):
    t = np.linspace(t0, t1, n)
    if curve_name == "line":
        
        P = np.array([params['x0'], params['y0'], params.get('z0', 0.0)], dtype=float)
        Q = np.array([params['x1'], params['y1'], params.get('z1', 0.0)], dtype=float)
        R = P[np.newaxis,:] + np.outer(t, (Q - P))
    elif curve_name == "circle":
        a = float(params.get('a', 1.0))
        R = np.column_stack((a*np.cos(t), a*np.sin(t), np.zeros_like(t)))
    elif curve_name == "ellipse":
        a = float(params.get('a', 2.0))
        b = float(params.get('b', 1.0))
        R = np.column_stack((a*np.cos(t), b*np.sin(t), np.zeros_like(t)))
    elif curve_name == "helix":
        a = float(params.get('a', 1.0))
        b = float(params.get('b', 0.2))
        R = np.column_stack((a*np.cos(t), a*np.sin(t), b*t))
    elif curve_name == "cycloid":
        a = float(params.get('a', 1.0))
        R = np.column_stack((a*(t - np.sin(t)), a*(1 - np.cos(t)), np.zeros_like(t)))
    elif curve_name == "twisted_cubic":
        R = np.column_stack((t, t**2, t**3))
    elif curve_name == "catenary":
        C = float(params.get('C', 1.0))
        R = np.column_stack((t, C*np.cosh(t/C), np.zeros_like(t)))
    elif curve_name == "hyperbola":
        R = np.column_stack((np.cosh(t), np.sinh(t), np.zeros_like(t)))
    elif curve_name == "tractrix":
        R = np.column_stack((t - np.tanh(t), 1/np.cosh(t), np.zeros_like(t)))
    else:
        raise ValueError("Unknown curve: " + curve_name)
    return t, R

def derivatives_numeric(R, t):
    dt = np.gradient(t)
    Rp = np.gradient(R, axis=0) / dt[:,np.newaxis]
    Rpp = np.gradient(Rp, axis=0) / dt[:,np.newaxis]
    Rppp = np.gradient(Rpp, axis=0) / dt[:,np.newaxis]
    return Rp, Rpp, Rppp

def curvature_torsion_from_R(R, t):
    Rp, Rpp, Rppp = derivatives_numeric(R, t)
    cross = np.cross(Rp, Rpp)
    cross_norm = np.linalg.norm(cross, axis=1)
    Rp_norm = np.linalg.norm(Rp, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        curvature = cross_norm / (Rp_norm**3)
    triple = np.einsum('ij,ij->i', cross, Rppp)
    denom = cross_norm**2
    with np.errstate(divide='ignore', invalid='ignore'):
        torsion = triple / denom
    curvature = [None if (not np.isfinite(x)) else float(x) for x in curvature]
    torsion = [None if (not np.isfinite(x)) else float(x) for x in torsion]
    return curvature, torsion

def arc_length_numeric(R, t):
    Rp = np.gradient(R, axis=0) / np.gradient(t)[:,np.newaxis]
    speed = np.linalg.norm(Rp, axis=1)
    dt = np.diff(t)
    if len(dt)==0:
        return [0.0]
    cum = np.zeros_like(speed)
    cum[1:] = np.cumsum((speed[:-1] + speed[1:]) * 0.5 * dt)
    return cum.tolist()

@app.route('/')
def index():
    return send_from_directory('static', 'index.html')

@app.route('/api/compute', methods=['POST'])
def api_compute():
    data = request.json
    curve = data.get('curve')
    params = data.get('params', {})
    t0 = float(data.get('t0', 0.0))
    t1 = float(data.get('t1', 2*math.pi))
    n = int(data.get('n', 400))
    try:
        t, R = numeric_curve_positions(curve, params, t0, t1, n)
    except Exception as e:
        return jsonify({'ok': False, 'error': str(e)}), 400

    curvature, torsion = curvature_torsion_from_R(R, t)
    arclen = arc_length_numeric(R, t)

    
    xs = R[:,0].tolist()
    ys = R[:,1].tolist()
    zs = R[:,2].tolist()

    
    sym_formula = symbolic_formula_for(curve, params)

    return jsonify({
        'ok': True,
        'curve': curve,
        't': t.tolist(),
        'x': xs, 'y': ys, 'z': zs,
        'curvature': curvature,
        'torsion': torsion,
        'arc_length': arclen,
        'symbolic': sym_formula
    })

def symbolic_formula_for(curve, params):
    t = sp.symbols('t')
    if curve == "line":
        P = sp.Matrix([params.get('x0',0), params.get('y0',0), params.get('z0',0)])
        Q = sp.Matrix([params.get('x1',1), params.get('y1',0), params.get('z1',0)])
        expr = P + t*(Q-P)
    elif curve == "circle":
        a = sp.symbols('a')
        expr = sp.Matrix([a*sp.cos(t), a*sp.sin(t)])
    elif curve == "ellipse":
        a,b = sp.symbols('a b')
        expr = sp.Matrix([a*sp.cos(t), b*sp.sin(t)])
    elif curve == "helix":
        a,b = sp.symbols('a b')
        expr = sp.Matrix([a*sp.cos(t), a*sp.sin(t), b*t])
    elif curve == "cycloid":
        a = sp.symbols('a')
        expr = sp.Matrix([a*(t-sp.sin(t)), a*(1-sp.cos(t))])
    elif curve == "twisted_cubic":
        expr = sp.Matrix([t, t**2, t**3])
    elif curve == "catenary":
        C = sp.symbols('C')
        expr = sp.Matrix([t, C*sp.cosh(t/C)])
    elif curve == "hyperbola":
        expr = sp.Matrix([sp.cosh(t), sp.sinh(t)])
    elif curve == "tractrix":
        expr = sp.Matrix([t - sp.tanh(t), sp.sech(t)])
    else:
        return ''
    return sp.latex(expr)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)
