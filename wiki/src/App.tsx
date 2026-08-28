import { useEffect, useState, type ReactNode } from 'react'
import { Link, NavLink, Route, Routes, useLocation } from 'react-router-dom'

const BRAIN_APP_URL = import.meta.env.VITE_BRAIN_APP_URL || 'http://127.0.0.1:8001'
const PUF_APP_URL = import.meta.env.VITE_PUF_APP_URL || 'http://127.0.0.1:8000'

const navItems = [
  ['Model', '/model'],
  ['Brain Delivery', '/brain-delivery'],
  ['Off-target Atlas', '/offtarget-atlas'],
  ['Engineering', '/engineering'],
  ['Software', '/software'],
  ['Resources', '/resources'],
]

function Brand() {
  return <Link className="brand" to="/" aria-label="REWIRE modeling home"><span className="brand-glyph" aria-hidden="true"><i /><i /><i /></span><span><b>REWIRE</b><small>MODELING</small></span></Link>
}

function Header() {
  const [open, setOpen] = useState(false)
  const location = useLocation()
  useEffect(() => setOpen(false), [location.pathname])
  return <header className="site-header"><Brand /><button className="menu-button" type="button" aria-expanded={open} aria-controls="main-nav" onClick={() => setOpen(!open)}>Menu <span aria-hidden="true">{open ? '×' : '≡'}</span></button><nav id="main-nav" className={open ? 'open' : ''} aria-label="Main navigation">{navItems.map(([label, path]) => <NavLink key={path} to={path}>{label}</NavLink>)}</nav></header>
}

function Footer() {
  return <footer className="site-footer"><div><Brand /><p>Two transparent models for safer, more testable APOE RNA-editing design.</p></div><div><b>Explore</b>{navItems.slice(0, 4).map(([label, path]) => <Link key={path} to={path}>{label}</Link>)}</div><div><b>Run</b><a href={BRAIN_APP_URL}>Brain Delivery App ↗</a><a href={PUF_APP_URL}>PUF Atlas App ↗</a><Link to="/software">Local setup</Link></div><div><b>Open science</b><p>Source code and content are prepared for the team repository and iGEM deployment.</p><span className="license">© 2026 · CC BY 4.0 content</span></div></footer>
}

function Layout({ children }: { children: ReactNode }) {
  const { pathname } = useLocation()
  useEffect(() => {
    window.scrollTo({ top: 0 })
  }, [pathname])
  return <><Header /><main>{children}</main><Footer /></>
}

function Eyebrow({ children, light = false }: { children: ReactNode; light?: boolean }) {
  return <p className={`eyebrow${light ? ' light' : ''}`}><span />{children}</p>
}

function ArrowLink({ to, children, external = false, className = '' }: { to: string; children: ReactNode; external?: boolean; className?: string }) {
  return external ? <a className={`arrow-link ${className}`} href={to}>{children}<b aria-hidden="true">↗</b></a> : <Link className={`arrow-link ${className}`} to={to}>{children}<b aria-hidden="true">→</b></Link>
}

function PageHero({ index, eyebrow, title, intro, tone = 'navy', children }: { index: string; eyebrow: string; title: ReactNode; intro: string; tone?: 'navy' | 'teal' | 'paper'; children?: ReactNode }) {
  return <section className={`page-hero tone-${tone}`}><div className="page-number">{index}</div><div><Eyebrow light={tone !== 'paper'}>{eyebrow}</Eyebrow><h1>{title}</h1><p>{intro}</p>{children}</div></section>
}

function ModelConstellation() {
  return <div className="constellation" aria-label="Connected Brain Delivery and PUF OffTarget models"><div className="model-node brain"><span>MODEL 01</span><strong>Brain<br />Delivery</strong><small>dose → BBB → editor</small></div><div className="bridge-line"><i /><b>candidate panel</b><i /></div><div className="model-node puf"><span>MODEL 02</span><strong>PUF<br />Atlas</strong><small>sequence → loci → priority</small></div><div className="question-node">Which design should<br />we test next?</div></div>
}

function Home() {
  return <>
    <section className="home-hero"><div className="hero-copy"><Eyebrow>REWIRE · iGEM modeling</Eyebrow><h1>Two models,<br />one engineering<br /><em>question.</em></h1><p>How can we move an APOE RNA editor toward the brain while identifying the transcriptome sites most worth validating?</p><div className="hero-actions"><ArrowLink to="/model" className="filled">Follow the model story</ArrowLink><ArrowLink to="/software">Run the software</ArrowLink></div></div><ModelConstellation /></section>
    <section className="question-band"><span>THE DESIGN QUESTION</span><blockquote>“Can a sequence-aware safety screen and a multiscale delivery model guide the same experimental decision?”</blockquote><p>Our answer is a transparent handoff: PUF-OffTarget Atlas produces a candidate panel; the Brain Delivery Model converts that evidence into an aggregate off-target prior.</p></section>
    <section className="model-cards"><article className="brain-card"><span className="card-index">01</span><Eyebrow>Systems model</Eyebrow><h2>Brain Delivery<br />Digital Twin</h2><p>A coupled ODE model links administration, systemic distribution, BBB transport, intracellular expression and APOE codon-specific editing.</p><ul><li>APOE3-like desired state</li><li>APOE2-like double-edit risk</li><li>Liver and blood exposure</li></ul><ArrowLink to="/brain-delivery">Explore the digital twin</ArrowLink></article><article className="puf-card"><span className="card-index">02</span><Eyebrow>Sequence model</Eyebrow><h2>PUF-OffTarget<br />Atlas</h2><p>A transcriptome-scale screen finds motif matches, maps genomic loci and ranks candidates using sequence, expression, RNA accessibility and consequence evidence.</p><ul><li>Human, mouse or custom transcriptomes</li><li>Binding-only and editor-fusion modes</li><li>Auditable candidate reports</li></ul><ArrowLink to="/offtarget-atlas">Explore the atlas</ArrowLink></article></section>
    <section className="principles"><div><Eyebrow>Our modeling stance</Eyebrow><h2>Useful because the uncertainty is visible.</h2></div><div className="principle-list"><article><b>01</b><h3>Mechanistic</h3><p>States and fluxes correspond to explicit biological hypotheses.</p></article><article><b>02</b><h3>Traceable</h3><p>Parameters, transformations and output files remain inspectable.</p></article><article><b>03</b><h3>Falsifiable</h3><p>Each model points to measurements that can confirm or revise it.</p></article><article><b>04</b><h3>Bounded</h3><p>Outputs are prioritization evidence and are not clinically calibrated.</p></article></div></section>
    <Disclaimer />
  </>
}

function WorkflowDiagram() {
  const steps = [['01', 'PUF sequence'], ['02', 'Transcriptome scan'], ['03', 'Candidate panel'], ['04', 'Delivery simulation'], ['05', 'Wet-lab choice']]
  return <ol className="workflow-diagram" aria-label="Cross-model workflow">{steps.map(([number, label], index) => <li key={number}><b>{number}</b><span>{label}</span>{index < steps.length - 1 && <i aria-hidden="true">→</i>}</li>)}</ol>
}

function ModelOverview() {
  return <><PageHero index="01" eyebrow="Integrated model" title={<>From sequence space<br />to system behavior.</>} intro="The two models operate at different scales. Their value comes from a deliberately narrow, inspectable interface rather than an opaque end-to-end prediction." />
    <section className="content-section intro-grid"><div><Eyebrow>Why two models?</Eyebrow><h2>Different questions need different abstractions.</h2></div><div><p className="lead">The Atlas asks where a PUF could bind. The digital twin asks how delivery, expression and editing evolve over time. Keeping them independent makes assumptions easier to test.</p><p>The shared CSV is a design artifact, not hidden server communication. A researcher can inspect, edit, archive and re-upload every candidate row.</p></div></section>
    <section className="dark-workflow"><Eyebrow light>Data contract</Eyebrow><h2>A small bridge across biological scales.</h2><WorkflowDiagram /><div className="contract-grid"><article><span>Atlas output</span><code>brain_candidate_panel.csv</code><p>Top 100 unique genomic loci, ordered by PUF prioritization score.</p></article><article><span>Model input</span><code>effective off-target pool</code><p>Accessibility-weighted site abundance, binding and catalytic context.</p></article><article><span>Decision output</span><code>benefit × risk × exposure</code><p>Comparable simulated outcomes across design, route and dose.</p></article></div></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Shared schema</Eyebrow><h2>Every transferred value has a named origin.</h2></div><div className="table-wrap"><table><thead><tr><th>Brain field</th><th>Atlas evidence</th><th>Interpretation</th></tr></thead><tbody><tr><td><code>initial_pool</code></td><td>expression score</td><td>relative available RNA abundance</td></tr><tr><td><code>binding_score</code></td><td>sequence score</td><td>relative PUF recognition prior</td></tr><tr><td><code>accessibility</code></td><td>RNAplfold score</td><td>relative motif exposure</td></tr><tr><td><code>context_score</code></td><td>consequence score</td><td>relative editing consequence prior</td></tr><tr><td><code>validation_priority</code></td><td>risk priority</td><td>human-readable triage label</td></tr></tbody></table></div><p className="callout">When expression or accessibility evidence is absent, the exporter uses a conservative value of 1.0 and records every fallback in metadata.</p></section>
    <Disclaimer />
  </>
}

function BrainPipeline() {
  const items = [['M1', 'Absorption', 'Adep · Alymph'], ['M2', 'Distribution', 'blood · liver · brain'], ['M3', 'BBB', 'EC · endosome · ISF'], ['M4', 'Expression', 'nucleus · mRNA · Pbrain'], ['M5', 'Editing', '112 · 158 · off-target']]
  return <div className="brain-pipeline">{items.map(([id, title, states], i) => <article key={id}><b>{id}</b><span>{title}</span><small>{states}</small>{i < items.length - 1 && <i>→</i>}</article>)}</div>
}

function BrainDelivery() {
  return <><PageHero index="02" eyebrow="Systems model" tone="navy" title={<>Brain Delivery<br />Digital Twin</>} intro="A sequence-aware, literature-informed ODE model connects administration to APOE RNA editing while keeping codon-specific benefit and risk separate."><ArrowLink to={BRAIN_APP_URL} external className="hero-launch">Launch Brain App</ArrowLink></PageHero>
    <section className="content-section"><div className="section-title"><Eyebrow>Five modules</Eyebrow><h2>One continuous state trajectory.</h2><p>Active brain editor abundance, <code>P_brain(t)</code>, is the interface between delivery and RNA editing.</p></div><BrainPipeline /></section>
    <section className="split-feature navy-block"><div><Eyebrow light>Biological correction</Eyebrow><h2>APOE codons 112 and 158 are not interchangeable.</h2><p>APOE4 begins Arg112/Arg158. Editing C112 alone yields an APOE3-like state; editing both C112 and C158 yields an APOE2-like state tracked as a risk proxy.</p></div><div className="isoform-grid"><article><span>START</span><b>Arg<sup>112</sup> / Arg<sup>158</sup></b><small>APOE4-like</small></article><article className="desired"><span>DESIRED</span><b>Cys<sup>112</sup> / Arg<sup>158</sup></b><small>APOE3-like</small></article><article className="risk"><span>RISK PROXY</span><b>Cys<sup>112</sup> / Cys<sup>158</sup></b><small>APOE2-like</small></article><article><span>MIXED</span><b>Arg<sup>112</sup> / Cys<sup>158</sup></b><small>non-target state</small></article></div></section>
    <section className="content-section equation-section"><div><Eyebrow>Core reaction</Eyebrow><h2>Binding, unbinding and catalysis stay explicit.</h2><p>Design features rescale kinetic priors. PUF score and accessibility affect association; UC context and editing-window distance affect catalysis.</p></div><div className="equation-card"><code>E + S ⇄ ES → E + P</code><dl><div><dt>k<sub>on</sub></dt><dd>PUF score · accessibility</dd></div><div><dt>k<sub>cat</sub></dt><dd>UC context · distance</dd></div><div><dt>off-target</dt><dd>PUF mismatch · background</dd></div></dl></div></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Readouts</Eyebrow><h2>A decision vector, not one magic score.</h2></div><div className="output-grid"><article><b>Benefit</b><span>APOE3-like fraction</span><p>C112-only editing trajectory.</p></article><article><b>Editing risk</b><span>APOE2-like proxy</span><p>Double-edit trajectory at 112 and 158.</p></article><article><b>Specificity</b><span>Off-target burden</span><p>Local bystander, PUF-mediated and deaminase background components.</p></article><article><b>Exposure</b><span>Liver AUC · blood Cmax</span><p>Systemic constraints around brain delivery.</p></article></div></section>
    <section className="action-band brain-action"><div><span>INTERACTIVE MODEL</span><h2>Compare design, route, dose and duration.</h2></div><ArrowLink to={BRAIN_APP_URL} external className="filled">Launch Brain Delivery App</ArrowLink></section><Disclaimer />
  </>
}

const bases = ['A', 'A', 'C', 'G', 'U', 'C', 'U', 'A', 'U', 'A']
function SequenceFigure() {
  return <div className="sequence-figure" aria-label="PUF motif with mismatch comparison"><div className="sequence-row query"><span>QUERY</span>{bases.map((base, i) => <b className={`base base-${base.toLowerCase()}`} key={`${base}-${i}`}>{base}</b>)}</div><div className="match-row"><span></span>{bases.map((_, i) => <i key={i}>{i === 6 ? '×' : '│'}</i>)}</div><div className="sequence-row hit"><span>HIT</span>{bases.map((base, i) => <b className={`base base-${(i === 6 ? 'c' : base).toLowerCase()}`} key={`${base}-${i}`}>{i === 6 ? 'C' : base}</b>)}</div><div className="sequence-stats"><span>1 mismatch</span><span>90% identity</span><span>candidate locus</span></div></div>
}

function OfftargetAtlas() {
  return <><PageHero index="03" eyebrow="Sequence model" tone="teal" title={<>PUF-OffTarget<br />Atlas</>} intro="A transcriptome-scale prioritization engine that turns a short PUF recognition sequence into an auditable list of possible off-target loci."><ArrowLink to={PUF_APP_URL} external className="hero-launch">Launch Atlas App</ArrowLink></PageHero>
    <section className="content-section atlas-intro"><div><Eyebrow>The question</Eyebrow><h2>Where else could this PUF bind?</h2><p className="lead">The Atlas searches 8–12 nt RNA motifs across annotated transcriptomes with 0–3 substitutions, then adds biological context that sequence identity alone cannot provide.</p></div><SequenceFigure /></section>
    <section className="score-stack"><div className="section-title"><Eyebrow light>Priority model</Eyebrow><h2>Four evidence layers, one ranked queue.</h2></div><div className="score-layers"><article><b>01</b><span>Sequence</span><strong>motif identity</strong><p>Mismatch count, position weights and optional substitution matrix.</p></article><article><b>02</b><span>Expression</span><strong>tissue evidence</strong><p>GTEx-style TPM summaries with stable Ensembl ID matching.</p></article><article><b>03</b><span>Structure</span><strong>RNA accessibility</strong><p>Optional RNAplfold unpaired probability around each motif.</p></article><article><b>04</b><span>Consequence</span><strong>annotation context</strong><p>CDS, UTR, splice proximity and potential editable bases.</p></article></div><p className="score-formula">priority = weighted evidence after renormalizing genuinely missing components</p></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Interpretation ladder</Eyebrow><h2>What a candidate does—and does not—mean.</h2></div><div className="evidence-ladder"><article><span>LEVEL 1</span><h3>Sequence candidate</h3><p>A transcript contains a sufficiently similar motif.</p></article><article><span>LEVEL 2</span><h3>Context-supported</h3><p>Expression, accessibility or annotation raises validation priority.</p></article><article><span>LEVEL 3</span><h3>Potential editable base</h3><p>An editor-compatible base falls inside the configured window.</p></article><article><span>REQUIRED</span><h3>Experimental evidence</h3><p>Binding or editing must still be measured.</p></article></div><p className="callout coral">A high Atlas score is a heuristic prioritization score—not a probability that binding or editing occurs.</p></section>
    <section className="action-band puf-action"><div><span>TRANSCRIPTOME SCREEN</span><h2>Start with a PUF sequence. Leave with a validation plan.</h2></div><ArrowLink to={PUF_APP_URL} external className="filled">Launch PUF-OffTarget Atlas</ArrowLink></section><Disclaimer />
  </>
}

function Engineering() {
  const cycles = [
    ['01', 'Structure the chain', 'Connected delivery, BBB, expression and editing states.'],
    ['02', 'Expose uncertainty', 'Recorded literature priors and ranges instead of false precision.'],
    ['03', 'Design experiments', 'Mapped uncertain parameters to qPCR, Western and amplicon readouts.'],
    ['04', 'Add sequence evidence', 'Replaced a hand-picked off-target pool with candidate sites.'],
    ['05', 'Fit transparently', 'Prepared an explainable fitting loop for sparse early data.'],
    ['06', 'Separate APOE states', 'Tracked desired C112-only editing and C112/C158 double editing.'],
  ]
  return <><PageHero index="04" eyebrow="Engineering" tone="paper" title={<>Design–Build–Test–Learn,<br />encoded in the model.</>} intro="The model changed because each iteration exposed a decision, a missing measurement or a biological simplification that mattered." />
    <section className="content-section"><div className="section-title"><Eyebrow>Six iterations</Eyebrow><h2>From a connected ODE to a sequence-aware digital twin.</h2></div><div className="cycle-list">{cycles.map(([id, title, text]) => <article key={id}><b>{id}</b><div><h3>{title}</h3><p>{text}</p></div></article>)}</div></section>
    <section className="dbtl-section"><div><Eyebrow light>Closed loop</Eyebrow><h2>Predictions become measurements; measurements revise priors.</h2></div><div className="dbtl-wheel" aria-label="Design Build Test Learn loop"><span className="design">DESIGN<small>choose construct</small></span><span className="build">BUILD<small>prepare variants</small></span><span className="test">TEST<small>collect readouts</small></span><span className="learn">LEARN<small>update parameters</small></span><i>↻</i></div></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Next measurements</Eyebrow><h2>The model tells us what information is valuable.</h2></div><div className="experiment-grid"><article><span>TIME COURSE</span><h3>24 / 48 / 72 h editing</h3><p>Fit target-editing kinetics and separate delivery delay from catalysis.</p></article><article><span>LINKED AMPLICON</span><h3>APOE112 + APOE158</h3><p>Resolve APOE3-like benefit from APOE2-like double-edit risk.</p></article><article><span>CONTROLS</span><h3>PUF-only · inactive · free editor</h3><p>Separate PUF mismatch risk from deaminase background.</p></article><article><span>EXPRESSION</span><h3>qPCR + Western</h3><p>Distinguish expression limitation from editing limitation.</p></article></div></section><Disclaimer />
  </>
}

function CodeBlock({ children }: { children: string }) { return <pre><code>{children}</code></pre> }

function Software() {
  return <><PageHero index="05" eyebrow="Software" title={<>Run the tools.<br />Inspect every handoff.</>} intro="The Wiki is static; both scientific tools run as independent local or hosted applications. No candidate data moves between servers automatically." />
    <section className="content-section app-showcase"><article className="app-card brain-app"><span>APP 01 · PORT 8001</span><h2>Brain Delivery Model</h2><p>Interactive single-run and dose-optimization workspace with local Plotly outputs.</p><CodeBlock>{`cd brain_delivery_model/XiaobxTim-iGEM-REWIRE-iterated-v4\npython -m pip install -r requirements.txt\npython run_brain_app.py`}</CodeBlock><ArrowLink to={BRAIN_APP_URL} external className="filled">Launch Brain App</ArrowLink></article><article className="app-card puf-app"><span>APP 02 · PORT 8000</span><h2>PUF-OffTarget Atlas</h2><p>Transcriptome registration, queued scans, candidate browser and downloadable reports.</p><CodeBlock>{`cd puf-offtarget-atlas\npython -m pip install -e ".[dev]"\npufscan web`}</CodeBlock><ArrowLink to={PUF_APP_URL} external className="filled">Launch Atlas App</ArrowLink></article></section>
    <section className="bridge-guide"><div><Eyebrow light>Cross-model workflow</Eyebrow><h2>Download. Inspect. Upload.</h2></div><ol><li><b>1</b><div><h3>Run an Atlas scan</h3><p>Choose a transcriptome and enter an 8–12 nt PUF sequence.</p></div></li><li><b>2</b><div><h3>Export the panel</h3><p>Download <code>brain_candidate_panel.csv</code> and its metadata audit file.</p></div></li><li><b>3</b><div><h3>Open Brain App</h3><p>Select a design, route and duration; attach the CSV in the optional bridge field.</p></div></li><li><b>4</b><div><h3>Compare outcomes</h3><p>Download JSON or CSV results and retain them with experimental records.</p></div></li></ol></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Guardrails</Eyebrow><h2>Designed for reproducible use.</h2></div><div className="guardrail-grid"><article><b>5 MB</b><span>panel upload limit</span></article><article><b>10,000</b><span>maximum candidate rows</span></article><article><b>0–1</b><span>validated score range</span></article><article><b>0</b><span>automatic server-to-server transfers</span></article></div></section><Disclaimer />
  </>
}

function Resources() {
  return <><PageHero index="06" eyebrow="Resources" tone="paper" title={<>Evidence and<br />reproducibility.</>} intro="The model is only as defensible as the trail from biological claim to parameter, source code, test and planned measurement." />
    <section className="content-section"><div className="section-title"><Eyebrow>Core evidence</Eyebrow><h2>References that shape model structure.</h2></div><div className="reference-list"><article><b>01</b><div><h3>CU-REWIRE architecture and editing window</h3><p>Han et al. (2022), <em>Nucleic Acids Research</em>. Supports PUF–APOBEC architecture, UC context, editing-window and 10-repeat specificity assumptions.</p><a href="https://doi.org/10.1093/nar/gkac713">DOI: 10.1093/nar/gkac713 ↗</a></div></article><article><b>02</b><div><h3>APOE isoform biology</h3><p>Codons 112 and 158 define the state logic used to distinguish APOE4-like, APOE3-like and APOE2-like RNA products.</p><a href="https://www.ncbi.nlm.nih.gov/gene/348">NCBI Gene: APOE ↗</a></div></article><article><b>03</b><div><h3>AAV systemic distribution and BBB concept</h3><p>Literature supports the module structure while route- and capsid-specific rates remain broad, uncalibrated priors.</p><a href="https://doi.org/10.1016/j.xphs.2023.10.005">Journal of Pharmaceutical Sciences ↗</a></div></article><article><b>04</b><div><h3>Transcript reference and annotation</h3><p>GENCODE transcript FASTA and comprehensive annotations anchor transcript-to-genome mapping.</p><a href="https://www.gencodegenes.org/">GENCODE ↗</a></div></article></div></section>
    <section className="repro-band"><div><Eyebrow light>Reproducibility map</Eyebrow><h2>Claim → code → test → artifact</h2></div><div className="repro-grid"><article><span>PARAMETERS</span><code>config/base_config.yaml</code><p>Model priors and optimization thresholds.</p></article><article><span>PROVENANCE</span><code>parameter_provenance.yaml</code><p>Source, rationale, uncertainty and next measurement.</p></article><article><span>TESTS</span><code>tests/</code><p>State logic, sequence mapping, panel bridge and web contracts.</p></article><article><span>OUTPUT</span><code>CSV · JSON · HTML</code><p>Portable results that remain inspectable outside each app.</p></article></div></section>
    <section className="content-section"><div className="section-title"><Eyebrow>Limitations</Eyebrow><h2>What the current models cannot establish.</h2></div><div className="limits"><p>They do not estimate a patient-specific therapeutic dose.</p><p>They do not prove PUF binding or RNA editing at a candidate locus.</p><p>They do not replace capsid-, species- and route-specific calibration.</p><p>They compress an off-target panel instead of simulating every RNA state.</p><p>They do not infer causal safety from a composite prioritization score.</p><p>They require linked wet-lab measurements to fit uncertain parameters.</p></div></section><Disclaimer />
  </>
}

function Disclaimer() {
  return <aside className="disclaimer"><span>MODEL BOUNDARY</span><p>All outputs are literature-informed hypotheses for design comparison. They are <strong>not clinically calibrated</strong>, not probabilities of biological events and not substitutes for experimental validation.</p></aside>
}

function NotFound() { return <section className="not-found"><Eyebrow>404</Eyebrow><h1>That coordinate is outside this atlas.</h1><ArrowLink to="/">Return home</ArrowLink></section> }

export default function App() {
  return <Layout><Routes><Route path="/" element={<Home />} /><Route path="/model" element={<ModelOverview />} /><Route path="/brain-delivery" element={<BrainDelivery />} /><Route path="/offtarget-atlas" element={<OfftargetAtlas />} /><Route path="/engineering" element={<Engineering />} /><Route path="/software" element={<Software />} /><Route path="/resources" element={<Resources />} /><Route path="*" element={<NotFound />} /></Routes></Layout>
}
