import { render, screen } from '@testing-library/react'
import { MemoryRouter } from 'react-router-dom'
import { describe, expect, it } from 'vitest'

import App from './App'

function renderRoute(route: string) {
  return render(
    <MemoryRouter initialEntries={[route]}>
      <App />
    </MemoryRouter>,
  )
}

describe('dual-model wiki routes', () => {
  it('introduces the connected modeling workflow', () => {
    renderRoute('/')
    expect(screen.getByRole('heading', { name: /two models, one engineering question/i })).toBeInTheDocument()
    expect(screen.getAllByText(/not clinically calibrated/i).length).toBeGreaterThan(0)
  })

  it.each([
    ['/model', /from sequence space to system behavior/i],
    ['/brain-delivery', /brain delivery digital twin/i],
    ['/offtarget-atlas', /puf-offtarget atlas/i],
    ['/engineering', /design–build–test–learn/i],
    ['/software', /run the tools/i],
    ['/resources', /evidence and reproducibility/i],
  ])('renders %s', (route, heading) => {
    renderRoute(route)
    expect(screen.getByRole('heading', { name: heading })).toBeInTheDocument()
  })

  it('uses configurable external app links', () => {
    renderRoute('/software')
    const links = screen.getAllByRole('link', { name: /launch/i })
    expect(links.some((link) => link.getAttribute('href')?.includes('8000'))).toBe(true)
    expect(links.some((link) => link.getAttribute('href')?.includes('8001'))).toBe(true)
  })
})
