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

document.addEventListener("DOMContentLoaded", () => {
  renderSequencePreview();
  updateEditorFields();
  bySelector("[data-sequence-input]")?.addEventListener("input", renderSequencePreview);
  document.querySelectorAll('input[name="mode"]').forEach((input) => input.addEventListener("change", updateEditorFields));
  document.addEventListener("click", (event) => {
    const row = event.target.closest("tr[data-href]");
    if (row && !event.target.closest("a,button,input,select")) window.location.assign(row.dataset.href);
  });
});
