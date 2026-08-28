import { readFile, stat } from 'node:fs/promises'
import { basename, join } from 'node:path'

const routes = ['', 'model', 'brain-delivery', 'offtarget-atlas', 'engineering', 'software', 'resources']
const htmlFiles = routes.map((route) => join('dist', route, 'index.html'))

for (const file of htmlFiles) {
  const info = await stat(file)
  if (!info.isFile()) throw new Error(`Missing static route entry: ${file}`)
}

const documents = await Promise.all(htmlFiles.map((file) => readFile(file, 'utf8')))
const index = documents[0]
const cssHref = index.match(/href="([^"]+\.css)"/)?.[1]
if (!cssHref) throw new Error('Production index does not reference a CSS bundle.')
const css = await readFile(join('dist', 'assets', basename(cssHref)), 'utf8')

const remoteAsset = /<(?:script|link)\b[^>]*(?:src|href)=["']https?:\/\//i
if (documents.some((document) => remoteAsset.test(document)) || /url\(\s*["']?https?:\/\//i.test(css)) {
  throw new Error('Remote runtime asset found in the production Wiki bundle.')
}

console.log(`audited ${htmlFiles.length} static routes; no remote runtime assets`)
