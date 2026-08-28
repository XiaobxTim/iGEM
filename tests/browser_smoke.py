"""Browser smoke test; run with the three local servers as documented in README."""

from playwright.sync_api import sync_playwright


def collect_errors(page: object, errors: list[str]) -> None:
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(str(error)))


with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True)
    errors: list[str] = []

    wiki = browser.new_page(viewport={"width": 1440, "height": 1000})
    collect_errors(wiki, errors)
    wiki.goto("http://127.0.0.1:5173/model/", wait_until="networkidle")
    assert wiki.get_by_role("heading", name="From sequence space to system behavior.").is_visible()

    mobile = browser.new_page(viewport={"width": 390, "height": 844})
    collect_errors(mobile, errors)
    mobile.goto("http://127.0.0.1:5173/offtarget-atlas/", wait_until="networkidle")
    assert mobile.get_by_role("heading", name="PUF-OffTarget Atlas").is_visible()
    mobile.get_by_role("button", name="Menu", exact=True).click()
    assert mobile.get_by_role("link", name="Software").is_visible()

    brain = browser.new_page(viewport={"width": 1440, "height": 1000})
    collect_errors(brain, errors)
    brain.goto("http://127.0.0.1:8001", wait_until="networkidle")
    brain.locator("#candidate-panel").set_input_files(
        {
            "name": "brain_candidate_panel.csv",
            "mimeType": "text/csv",
            "buffer": (
                b"site_id,gene,initial_pool,binding_score,accessibility,context_score,"
                b"validation_priority,notes\n"
                b"puf-e2e,GENE1,0.8,0.9,0.7,0.6,High priority,browser test\n"
            ),
        }
    )
    brain.select_option('select[name="duration"]', "24")
    brain.get_by_role("button", name="Run simulation →").click()
    brain.locator("#form-status").filter(has_text="Complete").wait_for(timeout=30_000)
    assert brain.get_by_text("Candidate panel applied:").is_visible()
    assert brain.locator("#primary-plot .plotly").count() == 1

    atlas = browser.new_page(viewport={"width": 1280, "height": 900})
    collect_errors(atlas, errors)
    atlas.goto("http://127.0.0.1:8000", wait_until="networkidle")
    assert atlas.get_by_role("heading", name="Trace every possible match.").is_visible()

    browser.close()
    if errors:
        raise AssertionError("Browser console errors:\n" + "\n".join(errors))

print("browser smoke test passed")
