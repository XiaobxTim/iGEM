import react from '@vitejs/plugin-react'
import { defineConfig } from 'vitest/config'

const teamSlug = process.env.VITE_TEAM_SLUG?.replace(/^\/+|\/+$/g, '')

export default defineConfig({
  base: teamSlug ? `/${teamSlug}/` : '/',
  plugins: [react()],
  test: {
    environment: 'jsdom',
    globals: true,
    setupFiles: './src/test/setup.ts',
  },
})
