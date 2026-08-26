const form = document.querySelector('#simulation-form');
const statusLine = document.querySelector('#form-status');
const runButton = document.querySelector('#run-button');
const results = document.querySelector('#results');
const doseInput = document.querySelector('input[name="dose"]');
const doseOutput = document.querySelector('#dose-output');
const panelInput = document.querySelector('#candidate-panel');
const fileLabel = document.querySelector('#file-label');
let latestResult = null;

const plotLayout = (yTitle) => ({
  margin: { l: 52, r: 18, t: 24, b: 48 },
  paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',
  font: { family: 'Inter, Arial, sans-serif', color: '#40566d', size: 11 },
  xaxis: { title: 'Time (h)', gridcolor: '#e4e8eb', zeroline: false },
  yaxis: { title: yTitle, gridcolor: '#e4e8eb', zeroline: false },
  legend: { orientation: 'h', y: 1.12 },
});

doseInput.addEventListener('input', () => { doseOutput.textContent = Number(doseInput.value).toFixed(1); });
panelInput.addEventListener('change', () => { fileLabel.textContent = panelInput.files[0]?.name || 'Choose candidate CSV'; });
document.querySelectorAll('input[name="mode"]').forEach((input) => input.addEventListener('change', () => {
  document.querySelector('#dose-control').style.opacity = input.value === 'optimize' && input.checked ? '.45' : '1';
}));

const fmt = (value, digits = 3) => Number(value).toLocaleString(undefined, { maximumFractionDigits: digits });
const metric = (label, value, note) => `<article class="metric"><span>${label}</span><strong>${value}</strong><small>${note}</small></article>`;

function renderSingle(data) {
  const m = data.metrics;
  document.querySelector('#result-title').textContent = 'Single-simulation readout';
  document.querySelector('#metric-grid').innerHTML = [
    metric('APOE3-like', `${fmt(100 * m.apoe3_like_fraction_final, 1)}%`, 'C112-only edited fraction'),
    metric('APOE2-like risk', `${fmt(100 * m.apoe2_like_fraction_final, 1)}%`, 'C112 + C158 proxy'),
    metric('Off-target burden', fmt(m.off_target_burden_final), 'effective final burden'),
    metric('Specificity index', fmt(m.specificity_index_final), 'on/off model ratio'),
    metric('Brain PUF peak', fmt(m.P_brain_peak), 'normalized abundance'),
    metric('Liver AUC', fmt(m.AUC_liver), 'normalized exposure·h'),
    metric('Blood Cmax', fmt(m.Cmax_blood), 'normalized exposure'),
    metric('Final on-target rate', fmt(m.on_target_editing_rate_final), 'aggregate editing rate'),
  ].join('');
  document.querySelector('#primary-plot-title').textContent = 'APOE editing state fractions';
  document.querySelector('#secondary-plot-title').textContent = 'Brain editor abundance';
  Plotly.react('primary-plot', [
    { x: data.series.time, y: data.series.apoe3_like, name: 'APOE3-like', mode: 'lines', line: { color: '#155eef', width: 3 } },
    { x: data.series.time, y: data.series.apoe2_like, name: 'APOE2-like risk', mode: 'lines', line: { color: '#f5a524', width: 3 } },
    { x: data.series.time, y: data.series.off_target, name: 'Off-target', mode: 'lines', line: { color: '#b42318', width: 2, dash: 'dot' } },
  ], plotLayout('Fraction / burden'), { responsive: true, displaylogo: false });
  Plotly.react('secondary-plot', [{ x: data.series.time, y: data.series.pbrain, name: 'Pbrain', fill: 'tozeroy', mode: 'lines', line: { color: '#155eef', width: 2 }, fillcolor: 'rgba(21,94,239,.14)' }], plotLayout('Pbrain'), { responsive: true, displaylogo: false });
}

function renderOptimization(data) {
  const rows = data.dose_scan;
  const best = data.minimum_feasible_dose;
  document.querySelector('#result-title').textContent = 'Dose-optimization landscape';
  document.querySelector('#metric-grid').innerHTML = [
    metric('Doses screened', rows.length, '0.1–10 normalized units'),
    metric('Feasible doses', rows.filter((row) => row.feasible).length, 'all configured constraints'),
    metric('Minimum feasible', best ? fmt(best.dose, 2) : 'None', best ? 'normalized dose units' : 'under current thresholds'),
    metric('Design', data.inputs.design_id, `${data.inputs.route} · ${data.inputs.duration_hours} h`),
  ].join('');
  document.querySelector('#primary-plot-title').textContent = 'Benefit and editing risk by dose';
  document.querySelector('#secondary-plot-title').textContent = 'Systemic exposure by dose';
  Plotly.react('primary-plot', [
    { x: rows.map(r => r.dose), y: rows.map(r => 100 * r.apoe3_like_fraction_final), name: 'APOE3-like %', mode: 'lines+markers', line: { color: '#155eef', width: 3 } },
    { x: rows.map(r => r.dose), y: rows.map(r => 100 * r.apoe2_like_fraction_final), name: 'APOE2-like risk %', mode: 'lines+markers', line: { color: '#f5a524', width: 3 } },
  ], { ...plotLayout('Final fraction (%)'), xaxis: { ...plotLayout('').xaxis, title: 'Normalized dose' } }, { responsive: true, displaylogo: false });
  Plotly.react('secondary-plot', [
    { x: rows.map(r => r.dose), y: rows.map(r => r.AUC_liver), name: 'Liver AUC', mode: 'lines+markers', line: { color: '#b42318' } },
    { x: rows.map(r => r.dose), y: rows.map(r => r.Cmax_blood), name: 'Blood Cmax', mode: 'lines+markers', line: { color: '#155eef' } },
  ], { ...plotLayout('Normalized exposure'), xaxis: { ...plotLayout('').xaxis, title: 'Normalized dose' } }, { responsive: true, displaylogo: false });
}

function renderPanelNote(data) {
  const panel = data.panel_summary;
  document.querySelector('#panel-note').innerHTML = panel.provided
    ? `<strong>Candidate panel applied:</strong> ${fmt(panel.n_sites, 0)} sites compressed into an effective pool of ${fmt(panel.effective_pool)}. The bridge preserves ranking evidence but does not model every RNA as a separate ODE state.`
    : '<strong>Base off-target prior used:</strong> no PUF candidate panel was uploaded. Run PUF-OffTarget Atlas for a transcriptome-informed comparison.';
}

form.addEventListener('submit', async (event) => {
  event.preventDefault(); runButton.disabled = true; statusLine.classList.remove('error'); statusLine.textContent = 'Solving coupled ODE modules…';
  try {
    const response = await fetch('/api/simulate', { method: 'POST', body: new FormData(form) });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || 'Simulation failed.');
    latestResult = data;
    if (data.mode === 'single') renderSingle(data); else renderOptimization(data);
    renderPanelNote(data); results.hidden = false; results.scrollIntoView({ behavior: 'smooth' });
    statusLine.textContent = 'Complete. Results remain in this browser unless downloaded.';
  } catch (error) { statusLine.textContent = error.message; statusLine.classList.add('error'); }
  finally { runButton.disabled = false; }
});

function download(name, content, type) {
  const link = document.createElement('a'); link.href = URL.createObjectURL(new Blob([content], { type })); link.download = name; link.click(); URL.revokeObjectURL(link.href);
}
document.querySelector('#download-json').addEventListener('click', () => { if (latestResult) download('brain_model_result.json', JSON.stringify(latestResult, null, 2), 'application/json'); });
document.querySelector('#download-csv').addEventListener('click', () => {
  if (!latestResult) return;
  const rows = latestResult.mode === 'optimize' ? latestResult.dose_scan : [latestResult.metrics];
  const keys = Object.keys(rows[0]); const csv = [keys.join(','), ...rows.map(row => keys.map(key => JSON.stringify(row[key] ?? '')).join(','))].join('\n');
  download('brain_model_result.csv', csv, 'text/csv');
});
