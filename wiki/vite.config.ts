import react from '@vitejs/plugin-react'
import { fileURLToPath } from 'node:url'
import { defineConfig } from 'vitest/config'

const teamSlug = process.env.VITE_TEAM_SLUG?.replace(/^\/+|\/+$/g, '')
const routes = ['', 'model', 'brain-delivery', 'offtarget-atlas', 'engineering', 'software', 'resources']

export default defineConfig({
  base: teamSlug ? `/${teamSlug}/` : '/',
  plugins: [react()],
  build: {
    rollupOptions: {
      input: Object.fromEntries(
        routes.map((route) => [
          route || 'home',
          fileURLToPath(new URL(`${route ? `${route}/` : ''}index.html`, import.meta.url)),
        ]),
      ),
    },
  },
  test: {
    environment: 'jsdom',
    globals: true,
    setupFiles: './src/test/setup.ts',
  },
})
