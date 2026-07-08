/// <reference types="vitest/config" />
import { defineConfig } from 'vitest/config'
import react from '@vitejs/plugin-react'

// Vite + React + TS. Tests run under vitest (node env; physics is pure math).
export default defineConfig({
  plugins: [react()],
  test: { environment: 'node', include: ['src/**/*.test.ts'] },
})
