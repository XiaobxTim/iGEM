const bySelector = (selector, root = document) => root.querySelector(selector);

function renderSequencePreview() {
  const input = bySelector("[data-sequence-input]");
  const preview = bySelector("[data-sequence-preview]");
  if (!input || !preview) return;
  const sequence = input.value.toUpperCase().replaceAll("T", "U").slice(0, 12);
  preview.replaceChildren(...Array.from(sequence, (letter) => {
    const base = document.createElement("span");
    base.className = `base base-${letter.toLowerCase()}`;
    base.textContent = letter;
    return base;
  }));
}

function updateEditorFields() {
  const selected = bySelector('input[name="mode"]:checked');
  document.querySelectorAll(".editor-only").forEach((field) => {
    field.hidden = selected?.value !== "editor_fusion";
  });
}

async function fetchPartial(url, targetSelector, outer = false) {
  const target = bySelector(targetSelector);
  if (!target) return;
  const response = await fetch(url, { headers: { "X-Requested-With": "fetch" } });
  if (!response.ok) throw new Error(`Partial request failed (${response.status})`);
  const html = await response.text();
  if (outer) target.outerHTML = html;
  else target.innerHTML = html;
}

function scheduleProgressPoll() {
  const progress = bySelector('#run-progress[hx-get]');
  if (!progress) return;
  window.setTimeout(async () => {
    try {
      await fetchPartial(progress.getAttribute('hx-get'), '#run-progress', true);
      scheduleProgressPoll();
    } catch (error) {
      console.error(error);
    }
  }, 2000);
}

function bindCandidateFilters() {
  const form = bySelector('form.filter-bar[hx-get]');
  if (!form) return;
  let timeout;
  const refresh = () => {
    window.clearTimeout(timeout);
    timeout = window.setTimeout(async () => {
      const query = new URLSearchParams(new FormData(form));
      try {
        await fetchPartial(`${form.getAttribute('hx-get')}?${query}`, '#candidate-table');
      } catch (error) {
        console.error(error);
      }
    }, 300);
  };
  form.addEventListener('change', refresh);
  form.addEventListener('input', refresh);
}

document.addEventListener("DOMContentLoaded", () => {
  renderSequencePreview();
  updateEditorFields();
  bySelector("[data-sequence-input]")?.addEventListener("input", renderSequencePreview);
  document.querySelectorAll('input[name="mode"]').forEach((input) => input.addEventListener("change", updateEditorFields));
  bindCandidateFilters();
  scheduleProgressPoll();
  document.addEventListener("click", (event) => {
    const partialButton = event.target.closest('button[hx-get][hx-target]');
    if (partialButton) {
      event.preventDefault();
      fetchPartial(
        partialButton.getAttribute('hx-get'),
        partialButton.getAttribute('hx-target'),
      ).catch((error) => console.error(error));
      return;
    }
    const row = event.target.closest("tr[data-href]");
    if (row && !event.target.closest("a,button,input,select")) window.location.assign(row.dataset.href);
  });
});
