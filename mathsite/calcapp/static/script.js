console.log("script.js loaded");

const modeSelect = document.getElementById("mode");
const curveSection = document.getElementById("curve-section");
const surfaceSection = document.getElementById("surface-section");

const quantitySelect = document.getElementById("quantity");

const curveSelect = document.getElementById("curve");
const surfaceSelect = document.getElementById("surface");
const curveParamsDiv = document.getElementById("curve-params");
const surfaceParamsDiv = document.getElementById("surface-params");

const plotDiv = document.getElementById("plot");
const infoDiv = document.getElementById("info");
const computeBtn = document.getElementById("compute");

/* NEW: tab elements */
const mainTabs = document.querySelectorAll(".main-tabs .tab");
const curveQuantityTabsContainer = document.getElementById("curve-quantity-tabs");
const surfaceQuantityTabsContainer = document.getElementById("surface-quantity-tabs");
const curveQuantityTabs = curveQuantityTabsContainer
  ? curveQuantityTabsContainer.querySelectorAll(".tab")
  : [];
const surfaceQuantityTabs = surfaceQuantityTabsContainer
  ? surfaceQuantityTabsContainer.querySelectorAll(".tab")
  : [];

/* Helper to set active class on tabs in a group */
function setActiveTab(tabs, activeTab) {
  tabs.forEach((tab) => {
    if (tab === activeTab) {
      tab.classList.add("active");
    } else {
      tab.classList.remove("active");
    }
  });
}

/* MAIN TAB HANDLERS (Curves / Surfaces) */
mainTabs.forEach((tab) => {
  tab.addEventListener("click", () => {
    const mode = tab.dataset.modeTab; // "curve" or "surface"

    // Visually mark the active main tab
    setActiveTab(mainTabs, tab);

    // Update hidden select + trigger existing change logic
    if (modeSelect) {
      modeSelect.value = mode;
      modeSelect.dispatchEvent(new Event("change"));
    }

    // Show relevant quantity subtabs
    if (mode === "curve") {
      if (curveQuantityTabsContainer) {
        curveQuantityTabsContainer.style.display = "flex";
      }
      if (surfaceQuantityTabsContainer) {
        surfaceQuantityTabsContainer.style.display = "none";
      }
    } else {
      if (curveQuantityTabsContainer) {
        curveQuantityTabsContainer.style.display = "none";
      }
      if (surfaceQuantityTabsContainer) {
        surfaceQuantityTabsContainer.style.display = "flex";
      }
    }
  });
});

/* QUANTITY SUBTAB HANDLER */
function handleQuantityTabClick(tab, groupTabs) {
  const quantity = tab.dataset.quantityTab || "";

  // Visually mark the active subtab in this group
  setActiveTab(groupTabs, tab);

  // Update hidden quantity select so backend logic stays the same
  if (quantitySelect) {
    quantitySelect.value = quantity;
  }

  // For surfaces, automatically toggle symbolic checkbox:
  // if a symbolic quantity is chosen, turn it on; for "Numeric only", off.
  const computeSymbolicCheckbox = document.getElementById("compute_symbolic");
  if (computeSymbolicCheckbox && modeSelect && modeSelect.value === "surface") {
    computeSymbolicCheckbox.checked = quantity !== "";
  }
}

/* Wire up curve quantity tabs */
curveQuantityTabs.forEach((tab) => {
  tab.addEventListener("click", () =>
    handleQuantityTabClick(tab, curveQuantityTabs)
  );
});

/* Wire up surface quantity tabs */
surfaceQuantityTabs.forEach((tab) => {
  tab.addEventListener("click", () =>
    handleQuantityTabClick(tab, surfaceQuantityTabs)
  );
});

/* Parameter templates */
const paramTemplates = {
  helix: [["a", "1.0"], ["b", "0.2"]],
  circle: [["a", "1.0"]],
  ellipse: [["a", "2.0"], ["b", "1.0"]],
  line: [
    ["x0", "0"],
    ["y0", "0"],
    ["z0", "0"],
    ["x1", "1"],
    ["y1", "0"],
    ["z1", "0"],
  ],
  cycloid: [["a", "1.0"]],
  twisted_cubic: [],
  catenary: [["C", "1.0"]],
  hyperbola: [],
  tractrix: [],
};

const surfaceTemplates = {
  sphere: [["r", "1.0"]],
  torus: [["R", "1.0"], ["r", "0.4"]],
  paraboloid: [["a", "1.0"]],
  custom_surface: [],
};

/* Curve parameter rendering */
function renderCurveParams(curve) {
  curveParamsDiv.innerHTML = "";

  // Custom curve: user defines x(t), y(t), z(t)
  if (curve === "custom_curve") {
    curveParamsDiv.innerHTML = `
      <label>Custom x(t), y(t), z(t)</label>
      <textarea id="curve_x" rows="2" placeholder="cos(t)"></textarea>
      <textarea id="curve_y" rows="2" placeholder="sin(t)"></textarea>
      <textarea id="curve_z" rows="2" placeholder="t"></textarea>
      <div style="font-size:12px;color:#666">Use Python-style: sin(t), cos(t), t**2, etc.</div>
    `;
    return;
  }

  // Special layout for line: x-range, y-range, z-range,
  // each on its own line, but with endpoints (0/1, etc.) inline.
  if (curve === "line") {
    const defaults = {};
    (paramTemplates.line || []).forEach(([name, val]) => {
      defaults[name] = val;
    });

    curveParamsDiv.innerHTML = `
      <div class="param-row">
        <label>x-range:</label>
        <input id="param_x0" value="${defaults.x0 ?? "0"}"> 
        <span class="inline-label">to</span>
        <input id="param_x1" value="${defaults.x1 ?? "1"}">
      </div>
      <div class="param-row">
        <label>y-range:</label>
        <input id="param_y0" value="${defaults.y0 ?? "0"}"> 
        <span class="inline-label">to</span>
        <input id="param_y1" value="${defaults.y1 ?? "0"}">
      </div>
      <div class="param-row">
        <label>z-range:</label>
        <input id="param_z0" value="${defaults.z0 ?? "0"}"> 
        <span class="inline-label">to</span>
        <input id="param_z1" value="${defaults.z1 ?? "0"}">
      </div>
    `;
    return;
  }

  // All other curves: one parameter per line
  const list = paramTemplates[curve] || [];
  if (!list.length) {
    curveParamsDiv.innerHTML = `<em style="color:#666;">No parameters (uses defaults).</em>`;
    return;
  }

  list.forEach(([name, val]) => {
    const row = document.createElement("div");
    row.className = "param-row";
    row.innerHTML = `
      <label>${name}</label>
      <input id="param_${name}" value="${val}">
    `;
    curveParamsDiv.appendChild(row);
  });
}

/* Surface parameter rendering */
function renderSurfaceParams(surface) {
  surfaceParamsDiv.innerHTML = "";
  if (surface === "custom_surface") {
    surfaceParamsDiv.innerHTML = `
      <label>Custom X(u,v), Y(u,v), Z(u,v)</label>
      <textarea id="surf_x" rows="2" placeholder="u*cos(v)"></textarea>
      <textarea id="surf_y" rows="2" placeholder="u*sin(v)"></textarea>
      <textarea id="surf_z" rows="2" placeholder="cos(u)+sin(v)"></textarea>
      <div style="font-size:12px;color:#666">Use u, v and Python-style functions.</div>
    `;
    return;
  }

  const list = surfaceTemplates[surface] || [];
  if (!list.length) {
    surfaceParamsDiv.innerHTML = `<em style="color:#666;">No parameters (uses defaults).</em>`;
    return;
  }

  list.forEach(([name, val]) => {
    const row = document.createElement("div");
    row.className = "param-row";
    row.innerHTML = `
      <label>${name}</label>
      <input id="surf_param_${name}" value="${val}">
    `;
    surfaceParamsDiv.appendChild(row);
  });
}

/* Mode switching (driven by hidden <select id="mode">) */
modeSelect.addEventListener("change", () => {
  const mode = modeSelect.value;
  if (mode === "curve") {
    curveSection.style.display = "block";
    surfaceSection.style.display = "none";
  } else {
    curveSection.style.display = "none";
    surfaceSection.style.display = "block";
  }
});

/* Curve/surface change listeners */
curveSelect.addEventListener("change", () => renderCurveParams(curveSelect.value));
surfaceSelect.addEventListener("change", () => renderSurfaceParams(surfaceSelect.value));

/* Initial render */
renderCurveParams(curveSelect.value);
renderSurfaceParams(surfaceSelect.value);

// Ensure initial hidden selects match default tabs
if (modeSelect) {
  modeSelect.value = "curve";
}
if (quantitySelect) {
  quantitySelect.value = "";
}

/* Compute handler */
computeBtn.addEventListener("click", async () => {
  computeBtn.disabled = true;
  computeBtn.textContent = "Computing...";
  infoDiv.textContent = "";

  const mode = modeSelect.value;

  try {
    if (mode === "curve") {
      const curve = curveSelect.value;
      let params = {};

      if (curve === "custom_curve") {
        params.exprs = {
          x: document.getElementById("curve_x").value.trim() || "0",
          y: document.getElementById("curve_y").value.trim() || "0",
          z: document.getElementById("curve_z").value.trim() || "0",
        };
      } else {
        const inputs = curveParamsDiv.querySelectorAll("input");
        inputs.forEach((inp) => {
          const name = inp.id.replace("param_", "");
          params[name] = inp.value;
        });
      }

      const t0 = parseFloat(document.getElementById("t0").value);
      const t1 = parseFloat(document.getElementById("t1").value);
      const n = parseInt(document.getElementById("n").value) || 400;

      const symbolic_quantity = quantitySelect.value; // curves: arc_length, reparam_arc_length, frenet

      const payload = {
        mode: "curve",
        curve,
        params,
        t0,
        t1,
        n,
        symbolic_quantity,
      };

      const resp = await fetch("/compute/", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      const text = await resp.text();
      let data;
      try {
        data = JSON.parse(text);
      } catch (e) {
        console.error("Failed to parse JSON:", text);
        throw new Error("Server did not return valid JSON.");
      }

      computeBtn.disabled = false;
      computeBtn.textContent = "Compute & Plot";

      if (!resp.ok || !data.ok) {
        infoDiv.textContent =
          "Error from server: " + (data && data.error ? data.error : resp.status);
        return;
      }

      const is3D = data.z.some((z) => Math.abs(z) > 1e-10);
      let trace;
      if (is3D) {
        trace = {
          type: "scatter3d",
          mode: "lines",
          x: data.x,
          y: data.y,
          z: data.z,
          line: { width: 4 },
        };
      } else {
        trace = {
          type: "scatter",
          mode: "lines",
          x: data.x,
          y: data.y,
          line: { width: 3 },
        };
      }

      Plotly.newPlot(plotDiv, [trace], {
        title: "Curve",
        margin: { l: 0, r: 0, t: 30, b: 0 },
        scene: { aspectmode: "auto" },
      });

      let html = "<h3>Curve Results</h3>";
      html += `<p>Points: ${data.x.length}</p>`;
      const idxs = [0, Math.floor(data.t.length / 2), data.t.length - 1];
      html +=
        "<table><tr><th>t</th><th>x</th><th>y</th><th>z</th><th>κ</th><th>τ</th><th>s</th></tr>";
      idxs.forEach((i) => {
        html += `<tr>
          <td>${data.t[i].toFixed(4)}</td>
          <td>${data.x[i].toFixed(4)}</td>
          <td>${data.y[i].toFixed(4)}</td>
          <td>${data.z[i].toFixed(4)}</td>
          <td>${data.curvature[i] === null ? "—" : data.curvature[i]}</td>
          <td>${data.torsion[i] === null ? "—" : data.torsion[i]}</td>
          <td>${data.arc_length[i].toFixed(4)}</td>
        </tr>`;
      });
      html += "</table>";

      if (data.symbolic && Object.keys(data.symbolic).length > 0) {
        html += "<h4>Symbolic results</h4>";
        for (const key in data.symbolic) {
          const latex = data.symbolic[key]; // already LaTeX from backend
          html += `<p><strong>${key}</strong>: <span class="math">\\[${latex}\\]</span></p>`;
        }
      }
      infoDiv.innerHTML = html;

      // Ask MathJax to typeset the new content
      if (window.MathJax && window.MathJax.typesetPromise) {
        MathJax.typesetPromise([infoDiv]).catch((err) => console.error(err));
      }
    } else {
      // Surface mode
      const surface = surfaceSelect.value;
      let params = {};

      if (surface === "custom_surface") {
        params.x = document.getElementById("surf_x").value.trim() || "0";
        params.y = document.getElementById("surf_y").value.trim() || "0";
        params.z = document.getElementById("surf_z").value.trim() || "0";
      } else {
        const inputs = surfaceParamsDiv.querySelectorAll("input");
        inputs.forEach((inp) => {
          const name = inp.id.replace("surf_param_", "");
          params[name] = inp.value;
        });
      }

      params.u0 = parseFloat(document.getElementById("u0").value);
      params.u1 = parseFloat(document.getElementById("u1").value);
      params.v0 = parseFloat(document.getElementById("v0").value);
      params.v1 = parseFloat(document.getElementById("v1").value);
      params.nu = parseInt(document.getElementById("nu").value) || 60;
      params.nv = parseInt(document.getElementById("nv").value) || 60;

      const compute_symbolic = document.getElementById("compute_symbolic").checked;
      const symbolic_quantity = quantitySelect.value; // surfaces: first_form, second_form, etc.

      const payload = {
        mode: "surface",
        surface,
        params,
        compute_symbolic,
        symbolic_quantity,
      };

      const resp = await fetch("/compute/", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      const text = await resp.text();
      let data;
      try {
        data = JSON.parse(text);
      } catch (e) {
        console.error("Failed to parse JSON:", text);
        throw new Error("Server did not return valid JSON.");
      }

      computeBtn.disabled = false;
      computeBtn.textContent = "Compute & Plot";

      if (!resp.ok || !data.ok) {
        infoDiv.textContent =
          "Error from server: " + (data && data.error ? data.error : resp.status);
        return;
      }

      const X = data.X;
      const Y = data.Y;
      const Z = data.Z;

      const flatX = [];
      const flatY = [];
      const flatZ = [];
      for (let i = 0; i < X.length; i++) {
        for (let j = 0; j < X[0].length; j++) {
          flatX.push(X[i][j]);
          flatY.push(Y[i][j]);
          flatZ.push(Z[i][j]);
        }
      }

      const mesh = [
        {
          type: "mesh3d",
          x: flatX,
          y: flatY,
          z: flatZ,
          opacity: 0.9,
        },
      ];

      Plotly.newPlot(plotDiv, mesh, {
        title: "Surface",
        margin: { l: 0, r: 0, t: 30, b: 0 },
        scene: { aspectmode: "auto" },
      });

      let html = "<h3>Surface Results</h3>";
      html += `<p>Resolution: ${params.nu} × ${params.nv}</p>`;

      if (data.symbolic && Object.keys(data.symbolic).length > 0) {
        html += "<h4>Symbolic results</h4>";
        for (const key in data.symbolic) {
          const latex = data.symbolic[key];
          html += `<p><strong>${key}</strong>: <span class="math">\\[${latex}\\]</span></p>`;
        }
      }
      infoDiv.innerHTML = html;

      // Ask MathJax to typeset the new content
      if (window.MathJax && window.MathJax.typesetPromise) {
        MathJax.typesetPromise([infoDiv]).catch((err) => console.error(err));
      }
    }
  } catch (err) {
    console.error(err);
    computeBtn.disabled = false;
    computeBtn.textContent = "Compute & Plot";
    infoDiv.textContent = "Error in browser: " + err.message;
  }
});
