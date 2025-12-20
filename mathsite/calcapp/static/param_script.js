/* Curve types */
const curveDefaults = {
    helix: { var: "t", x: "cos(t)", y: "sin(t)", z: "0.2*t" },
    circle: { var: "t", x: "cos(t)", y: "sin(t)", z: "0" },
    ellipse: { var: "t", x: "2*cos(t)", y: "sin(t)", z: "0" },
    line: { var: "t", x: "x0 + t*(x1-x0)", y: "y0 + t*(y1-y0)", z: "z0 + t*(z1-z0)" },
    cycloid: { var: "t", x: "t - sin(t)", y: "1 - cos(t)", z: "0" },
    twisted_cubic: { var: "t", x: "t", y: "t^2", z: "t^3" },
    catenary: { var: "t", x: "t", y: "cosh(t/2)", z: "0" },
    hyperbola: { var: "t", x: "cosh(t)", y: "sinh(t)", z: "0" },
    tractrix: { var: "t", x: "t - tanh(t)", y: "sech(t)", z: "0" }
};

/* Surface types */

const surfaceDefaults = {

    plane: {
        u: "u",
        v: "v",
        x: "u",
        y: "v",
        z: "u + v"
    },

    cylinder: {
        u: "u",
        v: "v",
        x: "cos(u)",
        y: "sin(u)",
        z: "v"
    },

    cone: {
        u: "u",
        v: "v",
        x: "v * cos(u)",
        y: "v * sin(u)",
        z: "v"
    },

    paraboloid: {
        u: "u",
        v: "v",
        x: "u",
        y: "v",
        z: "u^2 + v^2"
    },

    hyperbolic_paraboloid: {
        u: "u",
        v: "v",
        x: "u",
        y: "v",
        z: "u^2 - v^2"
    },

    sphere: {
        u: "u",
        v: "v",
        x: "cos(v) * sin(u)",
        y: "sin(u) * sin(v)",
        z: "cos(u)"
    },

    torus: {
        u: "u",
        v: "v",
        x: "(R + r * cos(v)) * cos(u)",
        y: "(R + r * cos(v)) * sin(u)",
        z: "r * sin(v)"
    },

    helicoid: {
        u: "u",
        v: "v",
        x: "v * cos(u)",
        y: "v * sin(u)",
        z: "u"
    },

    catenoid: {
        u: "u",
        v: "v",
        x: "cosh(v) * cos(u)",
        y: "cosh(v) * sin(u)",
        z: "v"
    },

    mobius: {
        u: "u",
        v: "v",
        x: "(1 + v * cos(u/2)) * cos(u)",
        y: "(1 + v * cos(u/2)) * sin(u)",
        z: "v * sin(u/2)"
    },

    klein: {
        u: "u",
        v: "v",
        x: "(cos(u) * (cos(u/2) * (sqrt(2)+cos(v)) + sin(u/2) * sin(v)))",
        y: "(sin(u) * (cos(u/2) * (sqrt(2)+cos(v)) + sin(u/2) * sin(v)))",
        z: "(sin(u/2) * (sqrt(2)+cos(v)) - cos(u/2) * sin(v))"
    },

    enneper: {
        u: "u",
        v: "v",
        x: "u - (u^3)/3 + u*v^2",
        y: "v - (v^3)/3 + v*u^2",
        z: "u^2 - v^2"
    }
};


/* Curve handler */
document.getElementById("curve").addEventListener("change", () => {
    const c = curve.value;

    if (c === "custom_curve") {
        var1.value = "";
        coord1.value = "";
        coord2.value = "";
        coord3.value = "";
        return;
    }

    const def = curveDefaults[c];
    var1.value = def.var;
    coord1.value = def.x;
    coord2.value = def.y;
    coord3.value = def.z;
});


/* Surface handler */
document.getElementById("surface").addEventListener("change", () => {
    const s = document.getElementById("surface").value;

    const usEl = document.getElementById("us");
    const vsEl = document.getElementById("vs");
    const usxEl = document.getElementById("usx");
    const usyEl = document.getElementById("usy");
    const uszEl = document.getElementById("usz");

    if (!usEl || !vsEl || !usxEl || !usyEl || !uszEl) {
        
        return;
    }

    if (s === "custom_surface") {
        usEl.value = "";
        vsEl.value = "";
        usxEl.value = "";
        usyEl.value = "";
        uszEl.value = "";
        return;
    }

    const def = surfaceDefaults[s];
    if (def) {
        usEl.value = def.u || "u";
        vsEl.value = def.v || "v";
        usxEl.value = def.x || "u";
        usyEl.value = def.y || "v";
        uszEl.value = def.z || "0";
    } else {
        
        usEl.value = "u";
        vsEl.value = "v";
        usxEl.value = "u";
        usyEl.value = "v";
        uszEl.value = "0";
    }
});


// Trigger defaults on load
curve.dispatchEvent(new Event("change"));
surface.dispatchEvent(new Event("change"));
