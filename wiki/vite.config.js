var _a;
import react from '@vitejs/plugin-react';
import { defineConfig } from 'vitest/config';
var teamSlug = (_a = process.env.VITE_TEAM_SLUG) === null || _a === void 0 ? void 0 : _a.replace(/^\/+|\/+$/g, '');
export default defineConfig({
    base: teamSlug ? "/".concat(teamSlug, "/") : '/',
    plugins: [react()],
    test: {
        environment: 'jsdom',
        globals: true,
        setupFiles: './src/test/setup.ts',
    },
});
