# REWIRE dual-model iGEM Wiki

English, model-focused static Wiki for the Brain Delivery Digital Twin and
PUF-OffTarget Atlas. The Wiki links to the two applications but does not embed
or exchange scientific data with them.

## Local development

```bash
cp .env.example .env.local
npm install
npm run dev
```

Open the URL printed by Vite (normally `http://127.0.0.1:5173`). The app URLs
can be changed in `.env.local`.

## Test and build

```bash
npm test
npm run lint
npm run build
npm run audit
```

The deployable static site is written to `dist/`. `dist/` and `node_modules/`
are intentionally ignored by Git.

## iGEM deployment

Set `VITE_TEAM_SLUG` to the path segment used by the team Wiki before building:

```bash
VITE_TEAM_SLUG=your-team-slug npm run build
```

Upload the contents of `dist/` through the team Wiki deployment workflow. Each
client route has its own static `index.html`, so direct links such as
`/model/` and `/software/` work without a catch-all server rewrite. The
bundle loads no remote fonts, scripts, stylesheets or images at runtime. The
external URLs on the Resources page are ordinary links and are not fetched by
the application.

For hosted model apps, also set `VITE_BRAIN_APP_URL` and `VITE_PUF_APP_URL` at
build time.
