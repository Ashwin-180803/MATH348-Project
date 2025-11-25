const curveSelect = document.getElementById('curve');
const paramsDiv = document.getElementById('params');
const computeBtn = document.getElementById('compute');
const plotDiv = document.getElementById('plot');
const infoDiv = document.getElementById('info');

const paramTemplates = {
  helix: [
    ['a','1.0'], ['b','0.2']
  ],
  circle: [['a','1.0']],
  ellipse: [['a','2.0'], ['b','1.0']],
  line: [['x0','0'],['y0','0'],['z0','0'],['x1','1'],['y1','1'],['z1','0']],
  cycloid: [['a','1.0']],
  twisted_cubic: [],
  catenary: [['C','1.0']],
  hyperbola: [],
  tractrix: []
};

function renderParams(curve) {
  paramsDiv.innerHTML = '';
  const list = paramTemplates[curve] || [];
  list.forEach(([k,v])=>{
    const row = document.createElement('div'); row.className='param-row';
    const label = document.createElement('label'); label.textContent = k;
    const input = document.createElement('input'); input.id = 'param_'+k; input.value = v;
    row.appendChild(label); row.appendChild(input); paramsDiv.appendChild(row);
  });
  if(list.length===0) {
    paramsDiv.innerHTML = '<em style="color:#666">No params — uses default values.</em>';
  }
}

function gatherParams(curve) {
  const list = paramTemplates[curve] || [];
  const params = {};
  list.forEach(([k])=>{
    const el = document.getElementById('param_'+k);
    if(el) params[k] = el.value;
  });
  return params;
}

curveSelect.addEventListener('change', ()=> renderParams(curveSelect.value));
renderParams(curveSelect.value);

computeBtn.addEventListener('click', async ()=>{
  computeBtn.disabled = true;
  computeBtn.textContent = 'Computing...';
  infoDiv.innerHTML = '';

  const curve = curveSelect.value;
  const params = gatherParams(curve);
  const t0 = parseFloat(document.getElementById('t0').value);
  const t1 = parseFloat(document.getElementById('t1').value);
  const n = parseInt(document.getElementById('n').value) || 400;

  const resp = await fetch('/api/compute', {
    method:'POST', headers: {'Content-Type':'application/json'},
    body: JSON.stringify({curve, params, t0, t1, n})
  });
  const j = await resp.json();
  computeBtn.disabled = false;
  computeBtn.textContent = 'Compute & Plot';
  if(!j.ok) {
    infoDiv.textContent = 'Error: ' + j.error;
    return;
  }

  // decide 2D or 3D
  const is3D = j.z.some(z => Math.abs(z) > 1e-12);
  const trace = {
    x: j.x, y: j.y, mode:'lines', name: j.curve
  };
  if(is3D) {
    trace.z = j.z;
  }

  const data = is3D ? [Object.assign({}, trace, {type:'scatter3d', line:{width:4}})] : [Object.assign({}, trace, {type:'scatter', line:{width:3}})];

  const layout = {
    title: j.curve + ' (interactive)',
    margin: {l:0,r:0,t:30,b:0},
    scene: {aspectmode:'auto'}
  };
  Plotly.newPlot(plotDiv, data, layout);

  // Show info: sample curvature/torsion/arc-length and LaTeX formula
  const nShow = Math.min(5, j.t.length);
  const sampleIdx = Array.from({length:nShow}, (_,i)=> Math.round(i*(j.t.length-1)/(nShow-1)));
  let html = '<h3>Results</h3>';
  html += '<p><strong>Symbolic parametric:</strong><br><code>\\(' + j.symbolic + '\\)</code></p>';
  html += '<p><strong>Sample values (first, mid, last):</strong></p>';
  html += '<table><tr><th>t</th><th>x</th><th>y</th><th>z</th><th>κ</th><th>τ</th><th>s</th></tr>';
  sampleIdx.forEach(i=>{
    html += `<tr>
      <td>${Number(j.t[i]).toFixed(4)}</td>
      <td>${Number(j.x[i]).toFixed(4)}</td>
      <td>${Number(j.y[i]).toFixed(4)}</td>
      <td>${Number(j.z[i]).toFixed(4)}</td>
      <td>${j.curvature[i]===null?'—':Number(j.curvature[i]).toExponential(3)}</td>
      <td>${j.torsion[i]===null?'—':Number(j.torsion[i]).toExponential(3)}</td>
      <td>${Number(j.arc_length[i]).toFixed(4)}</td>
    </tr>`;
  });
  html += '</table>';
  html += '<p>Tip: Drag & scroll inside the plot to rotate/zoom. Click & hold to pan.</p>';
  infoDiv.innerHTML = html;

  // render LaTeX using MathJax if available (we don't include it by default). Fallback: raw LaTeX shown.
});
