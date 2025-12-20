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


const curveParamDisplay = document.getElementById("curve-param-display");
const surfaceParamDisplay = document.getElementById("surface-param-display");


const mainTabs = document.querySelectorAll(".main-tabs .tab");
const curveQuantityTabsContainer = document.getElementById("curve-quantity-tabs");
const surfaceQuantityTabsContainer = document.getElementById("surface-quantity-tabs");
const curveQuantityTabs = curveQuantityTabsContainer
  ? curveQuantityTabsContainer.querySelectorAll(".tab")
  : [];
const surfaceQuantityTabs = surfaceQuantityTabsContainer
  ? surfaceQuantityTabsContainer.querySelectorAll(".tab")
  : [];


function setActiveTab(tabs, activeTab) {
  tabs.forEach((tab) => {
    if (tab === activeTab) {
      tab.classList.add("active");
    } else {
      tab.classList.remove("active");
    }
  });
}


mainTabs.forEach((tab) => {
  tab.addEventListener("click", () => {
    const mode = tab.dataset.modeTab; 

    
    setActiveTab(mainTabs, tab);

    
    if (modeSelect) {
      modeSelect.value = mode;
      modeSelect.dispatchEvent(new Event("change"));
    }

    
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


function handleQuantityTabClick(tab, groupTabs) {
  const quantity = tab.dataset.quantityTab || "";

  
  setActiveTab(groupTabs, tab);

  
  if (quantitySelect) {
    quantitySelect.value = quantity;
  }

  enforceCurveQuantityConstraints();

  const computeSymbolicCheckbox = document.getElementById("compute_symbolic");
  if (computeSymbolicCheckbox && modeSelect && modeSelect.value === "surface") {
    computeSymbolicCheckbox.checked = quantity !== "";
  }
}


curveQuantityTabs.forEach((tab) => {
  tab.addEventListener("click", () =>
    handleQuantityTabClick(tab, curveQuantityTabs)
  );
});


surfaceQuantityTabs.forEach((tab) => {
  tab.addEventListener("click", () =>
    handleQuantityTabClick(tab, surfaceQuantityTabs)
  );
});


function enforceCurveQuantityConstraints() {
  if (!curveSelect || !quantitySelect) return;

  const options = Array.from(curveSelect.options || []);
  const blockLine = quantitySelect.value === "frenet";
  const blockArcLength = quantitySelect.value === "reparam_arc_length";
  const arcLengthUnsupported = new Set([
    "ellipse",
    "tractrix",
    "cycloid",
    "twisted_cubic",
    "catenary",
    "hyperbola",
  ]); // lacks closed-form reparam or known issues

  options.forEach((opt) => {
    const shouldHide =
      (blockLine && opt.value === "line") ||
      (blockArcLength && arcLengthUnsupported.has(opt.value));
    opt.disabled = shouldHide;
    opt.hidden = shouldHide;
  });

  const selected = options.find((opt) => opt.value === curveSelect.value);
  if (selected && (selected.disabled || selected.hidden)) {
    const fallback = options.find((opt) => !opt.disabled && !opt.hidden);
    if (fallback) {
      curveSelect.value = fallback.value;
      curveSelect.dispatchEvent(new Event("change"));
    }
  }
}


const paramTemplates = {
  helix: [],
  circle: [],
  ellipse: [],
  line: [
    ["x0", "0"],
    ["y0", "0"],
    ["z0", "0"],
    ["x1", "1"],
    ["y1", "0"],
    ["z1", "0"],
  ],
  cycloid: [],
  twisted_cubic: [],
  catenary: [],
  hyperbola: [],
  tractrix: [],
};

function getCurveParamDefaults(curve) {
  const defaults = {};
  (paramTemplates[curve] || []).forEach(([name, val]) => {
    defaults[name] = val;
  });
  return defaults;
}

const surfaceTemplates = {
  sphere: [],
  torus: [],
  paraboloid: [],
  custom_surface: [],
};


function latexTypeset(target) {
  if (window.MathJax && window.MathJax.typesetPromise) {
    MathJax.typesetPromise([target]).catch((err) => console.error(err));
  }
}

function updateCurveParamDisplay() {
  if (!curveParamDisplay || !curveSelect) return;

  const curve = curveSelect.value;
  let latex = "";

  switch (curve) {
    case "helix":
      latex = "\\mathbf{X}(t) = \\big(a\\cos t,\\; a\\sin t,\\; bt\\big)";
      break;
    case "circle":
      latex = "\\mathbf{X}(t) = \\big(a\\cos t,\\; a\\sin t\\big)";
      break;
    case "ellipse":
      latex = "\\mathbf{X}(t) = \\big(a\\cos t,\\; b\\sin t\\big)";
      break;
    case "line":
      latex =
        "\\mathbf{X}(t) = \\big(x_0 + t(x_1 - x_0),\\; y_0 + t(y_1 - y_0),\\; z_0 + t(z_1 - z_0)\\big)";
      break;
    case "cycloid":
      latex =
        "\\mathbf{X}(t) = \\big(a(t - \\sin t),\\; a(1 - \\cos t)\\big)";
      break;
    case "twisted_cubic":
      latex = "\\mathbf{X}(t) = \\big(t,\\; t^{2},\\; t^{3}\\big)";
      break;
    case "catenary":
      latex = "\\mathbf{X}(t) = \\big(t,\\; C\\cosh(t/C)\\big)";
      break;
    case "hyperbola":
      latex = "\\mathbf{X}(t) = \\big(\\cosh t,\\; \\sinh t\\big)";
      break;
    case "tractrix":
      latex = "\\mathbf{X}(t) = \\big(t - \\tanh t,\\; \\operatorname{sech} t\\big)";
      break;
    case "custom_curve": {
      const coord1El = document.getElementById("coord1");
      const coord2El = document.getElementById("coord2");
      const coord3El = document.getElementById("coord3");
      const cx = coord1El?.value || "";
      const cy = coord2El?.value || "";
      const cz = coord3El?.value || "";

      if (!cx && !cy && !cz) {
        latex = "\\mathbf{X}(t) = \\big(\\,x(t),\\; y(t),\\; z(t)\\,\\big)";
      } else {
        const parts = [];
        parts.push(cx || "x(t)");
        parts.push(cy || "y(t)");
        if (cz !== "") parts.push(cz || "z(t)");
        latex = `\\mathbf{X}(t) = \\big(${parts.join(",\\; ")}\\big)`;
      }
      break;
    }
    default:
      latex = "";
  }

  if (latex) {
    curveParamDisplay.innerHTML = `<span class="math">\\(${latex}\\)</span>`;
  } else {
    curveParamDisplay.innerHTML = "";
  }

  latexTypeset(curveParamDisplay);
}

function updateSurfaceParamDisplay() {
  if (!surfaceSelect) return;

  
  const varU = document.getElementById("us")?.value?.trim() || "u";
  const varV = document.getElementById("vs")?.value?.trim() || "v";

  
  const usxInput = document.getElementById("usx");
  const usyInput = document.getElementById("usy");
  const uszInput = document.getElementById("usz");

  if (usxInput) usxInput.placeholder = `X(${varU},${varV})`;
  if (usyInput) usyInput.placeholder = `Y(${varU},${varV})`;
  if (uszInput) uszInput.placeholder = `Z(${varU},${varV})`;

  
}


function renderCurveParams(curve) {
  curveParamsDiv.innerHTML = "";

  
  if (curve === "custom_curve") {
    curveParamsDiv.innerHTML = `<em style="color:#666;">Uses Defaults, unless changed.</em>`;
    updateCurveParamDisplay();
    return;
  }


  if (curve === "line") {
    curveParamsDiv.innerHTML = `<em style="color:#666;">Uses fixed defaults.</em>`;
    updateCurveParamDisplay();
    return;
  }

  
  const list = paramTemplates[curve] || [];
  if (!list.length) {
    curveParamsDiv.innerHTML = `<em style="color:#666;">Uses Defaults, unless changed.</em>`;
    updateCurveParamDisplay();
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

  updateCurveParamDisplay();
}


function renderSurfaceParams(surface) {
  surfaceParamsDiv.innerHTML = "";
  if (surface === "custom_surface") {
    
    const us = document.getElementById("us");
    const vs = document.getElementById("vs");
    const usx = document.getElementById("usx");
    const usy = document.getElementById("usy");
    const usz = document.getElementById("usz");
    [us, vs, usx, usy, usz].forEach((el) => {
      if (el) el.addEventListener("input", updateSurfaceParamDisplay);
    });

    surfaceParamsDiv.innerHTML = `<em style="color:#666;">Uses Defaults, unless changed.</em>`;
    updateSurfaceParamDisplay();
    return;
  }

  const list = surfaceTemplates[surface] || [];
  if (!list.length) {
    
    const us = document.getElementById("us");
    const vs = document.getElementById("vs");
    const usx = document.getElementById("usx");
    const usy = document.getElementById("usy");
    const usz = document.getElementById("usz");
    [us, vs, usx, usy, usz].forEach((el) => {
      if (el) el.addEventListener("input", updateSurfaceParamDisplay);
    });

    surfaceParamsDiv.innerHTML = `<em style="color:#666;">Uses Defaults, unless changed.</em>`;
    updateSurfaceParamDisplay();
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

  
  const us = document.getElementById("us");
  const vs = document.getElementById("vs");
  const usx = document.getElementById("usx");
  const usy = document.getElementById("usy");
  const usz = document.getElementById("usz");
  [us, vs, usx, usy, usz].forEach((el) => {
    if (el) el.addEventListener("input", updateSurfaceParamDisplay);
  });

  updateSurfaceParamDisplay();
}


modeSelect.addEventListener("change", () => {
  const mode = modeSelect.value;
  if (mode === "curve") {
    curveSection.style.display = "block";
    surfaceSection.style.display = "none";
  } else {
    curveSection.style.display = "none";
    surfaceSection.style.display = "block";
    
    if (surfaceSelect) {
      setTimeout(() => {
        surfaceSelect.dispatchEvent(new Event("change"));
      }, 10);
    }
  }
});


if (curveSelect) {
  curveSelect.addEventListener("change", () => {
    renderCurveParams(curveSelect.value);
    updateCurveParamDisplay();
  });
}
if (surfaceSelect) {
  surfaceSelect.addEventListener("change", () => {
    renderSurfaceParams(surfaceSelect.value);
    updateSurfaceParamDisplay();
  });
}


if (curveSelect) renderCurveParams(curveSelect.value);
if (surfaceSelect) renderSurfaceParams(surfaceSelect.value);
updateCurveParamDisplay();
updateSurfaceParamDisplay();
enforceCurveQuantityConstraints();


if (modeSelect) {
  modeSelect.value = "curve";
}
if (quantitySelect) {
  quantitySelect.value = "";
}
if (quantitySelect) {
  quantitySelect.addEventListener("change", enforceCurveQuantityConstraints);
}


if (computeBtn) {
  computeBtn.addEventListener("click", async () => {
    if (!modeSelect) return;
    computeBtn.disabled = true;
    computeBtn.textContent = "Computing...";
    if (infoDiv) infoDiv.textContent = "";

    const mode = modeSelect.value;

  try {
    if (mode === "curve") {
      const curve = curveSelect.value;
      let params = {};

      
      const var1El = document.getElementById("var1");
      const coord1El = document.getElementById("coord1");
      const coord2El = document.getElementById("coord2");
      const coord3El = document.getElementById("coord3");
      
      params.var = var1El && var1El.value.trim() ? var1El.value.trim() : "t";
      params.exprs = {
        x: coord1El && coord1El.value.trim() ? coord1El.value.trim() : "0",
        y: coord2El && coord2El.value.trim() ? coord2El.value.trim() : "0",
        z: coord3El && coord3El.value.trim() ? coord3El.value.trim() : "0",
      };

      if (curve === "line") {
        params = { ...params, ...getCurveParamDefaults("line") };
      } else if (curve === "custom_curve") {
        const inputs = curveParamsDiv.querySelectorAll("input");
        inputs.forEach((inp) => {
          const name = inp.id.replace("param_", "");
          if (name && inp.value) {
            params[name] = inp.value;
          }
        });
      }

      const t0El = document.getElementById("t0");
      const t1El = document.getElementById("t1");
      const nEl = document.getElementById("n");
      const t0 = t0El ? t0El.value.trim() : "0";
      const t1 = t1El ? t1El.value.trim() : "2*pi";
      const n = nEl ? parseInt(nEl.value) || 400 : 400;

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
        if (infoDiv) {
          infoDiv.textContent = "Error: Server did not return valid JSON. Please try again.";
        }
        if (computeBtn) {
          computeBtn.disabled = false;
          computeBtn.textContent = "Compute & Plot";
        }
        return;
      }

      if (computeBtn) {
        computeBtn.disabled = false;
        computeBtn.textContent = "Compute & Plot";
      }

      if (!resp.ok || !data || !data.ok) {
        const errorMsg = (data && data.error) ? data.error : `Server error (status: ${resp.status})`;
        if (infoDiv) {
          infoDiv.textContent = "Error from server: " + errorMsg;
        }
        return;
      }

      
      if (!data.x || !data.y || !data.z || !data.t ||
          !Array.isArray(data.x) || !Array.isArray(data.y) ||
          !Array.isArray(data.z) || !Array.isArray(data.t) ||
          data.x.length === 0 || data.y.length === 0 ||
          data.z.length === 0 || data.t.length === 0) {
        if (infoDiv) {
          infoDiv.textContent = "Error: Invalid response from server - missing or invalid data arrays";
        }
        return;
      }

      const is3D = data.z.some((z) => Math.abs(z) > 1e-10);
      let trace;
      let layout;
      
      if (is3D) {
        trace = {
          type: "scatter3d",
          mode: "lines",
          x: data.x,
          y: data.y,
          z: data.z,
          line: { width: 4, color: "#4f46e5" },
        };
        
        
        const xRange = [Math.min(...data.x), Math.max(...data.x)];
        const yRange = [Math.min(...data.y), Math.max(...data.y)];
        const zRange = [Math.min(...data.z), Math.max(...data.z)];
        
        const xPadding = (xRange[1] - xRange[0]) * 0.1 || 1;
        const yPadding = (yRange[1] - yRange[0]) * 0.1 || 1;
        const zPadding = (zRange[1] - zRange[0]) * 0.1 || 1;
        
        layout = {
          title: "Curve",
          margin: { l: 0, r: 0, t: 30, b: 0 },
          scene: {
            aspectmode: "auto",
            xaxis: { range: [xRange[0] - xPadding, xRange[1] + xPadding] },
            yaxis: { range: [yRange[0] - yPadding, yRange[1] + yPadding] },
            zaxis: { range: [zRange[0] - zPadding, zRange[1] + zPadding] },
            camera: {
              eye: { x: 1.5, y: 1.5, z: 1.5 }
            }
          },
        };
      } else {
        trace = {
          type: "scatter",
          mode: "lines",
          x: data.x,
          y: data.y,
          line: { width: 3, color: "#4f46e5" },
        };
        
        
        const xRange = [Math.min(...data.x), Math.max(...data.x)];
        const yRange = [Math.min(...data.y), Math.max(...data.y)];
        
        const xPadding = (xRange[1] - xRange[0]) * 0.1 || 1;
        const yPadding = (yRange[1] - yRange[0]) * 0.1 || 1;
        
        layout = {
          title: "Curve",
          margin: { l: 50, r: 20, t: 30, b: 50 },
          xaxis: {
            range: [xRange[0] - xPadding, xRange[1] + xPadding],
            title: "x",
            showgrid: true,
            zeroline: true,
          },
          yaxis: {
            range: [yRange[0] - yPadding, yRange[1] + yPadding],
            title: "y",
            showgrid: true,
            zeroline: true,
            scaleanchor: "x",
            scaleratio: 1,
          },
        };
      }

      Plotly.newPlot(plotDiv, [trace], layout);

      let html = "<h3>Curve Results</h3>";
      html += `<p>Points: ${data.x.length}</p>`;
      const idxs = [0, Math.floor(data.t.length / 2), data.t.length - 1];
      html +=
        "<table><tr><th>t</th><th>x</th><th>y</th><th>z</th><th>κ</th><th>τ</th><th>s</th></tr>";
      idxs.forEach((i) => {
        const kappa = data.curvature && data.curvature[i] !== null ? data.curvature[i].toFixed(4) : "—";
        const tau = data.torsion && data.torsion[i] !== null ? data.torsion[i].toFixed(4) : "—";
        const arc = data.arc_length && data.arc_length[i] !== null ? data.arc_length[i].toFixed(4) : "—";
        html += `<tr>
          <td>${data.t[i].toFixed(4)}</td>
          <td>${data.x[i].toFixed(4)}</td>
          <td>${data.y[i].toFixed(4)}</td>
          <td>${data.z[i].toFixed(4)}</td>
          <td>${kappa}</td>
          <td>${tau}</td>
          <td>${arc}</td>
        </tr>`;
      });
      html += "</table>";

      if (data.symbolic && Object.keys(data.symbolic).length > 0) {
        html += "<h4>Symbolic results</h4>";
        
        // Display computation steps if available
        if (data.symbolic.computation_steps && Array.isArray(data.symbolic.computation_steps)) {
          html += "<div style='margin-bottom: 15px; padding: 10px; background-color: #f3f4f6; border-left: 4px solid #3b82f6;'>";
          html += "<strong>Computation Steps:</strong><ol style='margin: 10px 0; padding-left: 20px;'>";
          data.symbolic.computation_steps.forEach(step => {
            html += `<li style='margin: 5px 0;'>${step}</li>`;
          });
          html += "</ol></div>";
        }
        
        for (const key in data.symbolic) {
          // Skip computation_steps as we already displayed it
          if (key === "computation_steps") continue;
          
          const value = data.symbolic[key]; // LaTeX from backend or error message
          
          // Check if this is an error field (should be displayed as plain text)
          if (key.endsWith("_error") || key.endsWith("_warning") || key.endsWith("_note")) {
            html += `<p><strong>${key.replace(/_/g, " ")}</strong>: <span style="color: #dc2626; font-style: italic;">${value}</span></p>`;
          } else {
            // Regular LaTeX content
            html += `<p><strong>${key}</strong>: <span class="math">\\[${value}\\]</span></p>`;
          }
        }
      }
      infoDiv.innerHTML = html;

      
      latexTypeset(infoDiv);
    } else {
      // Surface mode
      const surface = surfaceSelect.value;
      let params = {};

      
      const usEl = document.getElementById("us");
      const vsEl = document.getElementById("vs");
      const usxEl = document.getElementById("usx");
      const usyEl = document.getElementById("usy");
      const uszEl = document.getElementById("usz");
      params.u = usEl ? usEl.value.trim() || "u" : "u";
      params.v = vsEl ? vsEl.value.trim() || "v" : "v";
      params.x = usxEl ? usxEl.value.trim() : "0";
      params.y = usyEl ? usyEl.value.trim() : "0";
      params.z = uszEl ? uszEl.value.trim() : "0";

      
      const inputs = surfaceParamsDiv.querySelectorAll("input");
      inputs.forEach((inp) => {
        const name = inp.id.replace("surf_param_", "");
        params[name] = inp.value;
      });

      const u0El = document.getElementById("u0");
      const u1El = document.getElementById("u1");
      const v0El = document.getElementById("v0");
      const v1El = document.getElementById("v1");
      const nuEl = document.getElementById("nu");
      const nvEl = document.getElementById("nv");
      params.u0 = u0El ? u0El.value.trim() : "0";
      params.u1 = u1El ? u1El.value.trim() : "2*pi";
      params.v0 = v0El ? v0El.value.trim() : "0";
      params.v1 = v1El ? v1El.value.trim() : "2*pi";
      params.nu = nuEl ? parseInt(nuEl.value) || 60 : 60;
      params.nv = nvEl ? parseInt(nvEl.value) || 60 : 60;

      const compute_symbolic = document.getElementById("compute_symbolic").checked;
      const symbolic_quantity = quantitySelect.value; // surfaces: first_form, second_form, gaussian_curvature, mean_curvature, principal_curvatures, ...

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
        if (infoDiv) {
          infoDiv.textContent = "Error: Server did not return valid JSON. Please try again.";
        }
        if (computeBtn) {
          computeBtn.disabled = false;
          computeBtn.textContent = "Compute & Plot";
        }
        return;
      }

      if (computeBtn) {
        computeBtn.disabled = false;
        computeBtn.textContent = "Compute & Plot";
      }

      if (!resp.ok || !data || !data.ok) {
        const errorMsg = (data && data.error) ? data.error : `Server error (status: ${resp.status})`;
        if (infoDiv) {
          infoDiv.textContent = "Error from server: " + errorMsg;
        }
        return;
      }

      
      if (!data.X || !data.Y || !data.Z ||
          !Array.isArray(data.X) || !Array.isArray(data.Y) || !Array.isArray(data.Z) ||
          data.X.length === 0 || data.Y.length === 0 || data.Z.length === 0) {
        if (infoDiv) {
          infoDiv.textContent = "Error: Invalid response from server - missing or invalid surface data arrays";
        }
        return;
      }

      const X = data.X;
      const Y = data.Y;
      const Z = data.Z;

      
      const quantity = symbolic_quantity;
      let colorField = Z;
      let colorbarTitle = "Height (z)";

      if (quantity === "gaussian_curvature" && data.K) {
        colorField = data.K;
        colorbarTitle = "Gaussian curvature K";
      } else if (quantity === "mean_curvature" && data.H) {
        colorField = data.H;
        colorbarTitle = "Mean curvature H";
      } else if (quantity === "principal_curvatures" && data.k1) {
        
        colorField = data.k1;
        colorbarTitle = "Principal curvature k₁";
      }

      const surfaceTrace = {
        type: "surface",
        x: X,
        y: Y,
        z: Z,
        surfacecolor: colorField,
        colorscale: "Viridis",
        showscale: true,
        colorbar: {
          title: colorbarTitle,
          thickness: 15,
        },
      };

      Plotly.newPlot(plotDiv, [surfaceTrace], {
        title: "Surface",
        margin: { l: 0, r: 0, t: 30, b: 0 },
        scene: { aspectmode: "auto" },
      });

      let html = "<h3>Surface Results</h3>";
      html += `<p>Resolution: ${params.nu} × ${params.nv}</p>`;

      if (data.symbolic && Object.keys(data.symbolic).length > 0) {
        html += "<h4>Symbolic results</h4>";
        
        // Display computation steps if available
        if (data.symbolic.computation_steps && Array.isArray(data.symbolic.computation_steps)) {
          html += "<div style='margin-bottom: 15px; padding: 10px; background-color: #f3f4f6; border-left: 4px solid #3b82f6;'>";
          html += "<strong>Computation Steps:</strong><ol style='margin: 10px 0; padding-left: 20px;'>";
          data.symbolic.computation_steps.forEach(step => {
            html += `<li style='margin: 5px 0;'>${step}</li>`;
          });
          html += "</ol></div>";
        }
        
        for (const key in data.symbolic) {
          // Skip computation_steps as we already displayed it
          if (key === "computation_steps") continue;
          
          const value = data.symbolic[key]; // LaTeX from backend or error message
          
          // Check if this is an error field (should be displayed as plain text)
          if (key.endsWith("_error") || key.endsWith("_warning") || key.endsWith("_note")) {
            html += `<p><strong>${key.replace(/_/g, " ")}</strong>: <span style="color: #dc2626; font-style: italic;">${value}</span></p>`;
          } else {
            // Regular LaTeX content
            html += `<p><strong>${key}</strong>: <span class="math">\\[${value}\\]</span></p>`;
          }
        }
      }
      infoDiv.innerHTML = html;

      latexTypeset(infoDiv);
    }
  } catch (err) {
    console.error(err);
    if (computeBtn) {
      computeBtn.disabled = false;
      computeBtn.textContent = "Compute & Plot";
    }
    if (infoDiv) infoDiv.textContent = "Error in browser: " + err.message;
  }
  });
}