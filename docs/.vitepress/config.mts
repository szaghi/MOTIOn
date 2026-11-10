import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'

export default withMermaid({
  title: 'MOTIOn Documentation',
  base: '/MOTIOn/',
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/MOTIOn/favicon.svg', media: '(prefers-color-scheme: light)' }],
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/MOTIOn/favicon-dark.svg', media: '(prefers-color-scheme: dark)' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '32x32', href: '/MOTIOn/favicon-32.png' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '16x16', href: '/MOTIOn/favicon-16.png' }],
    ['link', { rel: 'apple-touch-icon', sizes: '180x180', href: '/MOTIOn/apple-touch-icon.png' }],
    ['meta', { name: 'theme-color', content: '#eff6f4', media: '(prefers-color-scheme: light)' }],
    ['meta', { name: 'theme-color', content: '#041210', media: '(prefers-color-scheme: dark)' }],
  ],
  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      'fortran': 'fortran-free-form',
      'f90': 'fortran-free-form',
      'f95': 'fortran-free-form',
      'f03': 'fortran-free-form',
      'f08': 'fortran-free-form',
      'f77': 'fortran-fixed-form',
    },
  },
  themeConfig: {
    logo: { light: '/logo-nav.svg', dark: '/logo-nav-dark.svg', alt: 'MOTIOn' },
    siteTitle: 'MOTIOn',
    nav: [
      { text: 'Home', link: '/' },
      {
        text: 'Guide',
        items: [
          { text: 'About',             link: '/guide/' },
          { text: 'Features',          link: '/guide/features' },
          { text: 'Installation',      link: '/guide/installation' },
          { text: 'Usage',             link: '/guide/usage' },
          { text: 'API Reference',     link: '/guide/api-reference' },
          { text: 'Contributing',      link: '/guide/contributing' },
          { text: 'Coverage Analysis', link: '/guide/coverage-analysis' },
          { text: 'Changelog',         link: '/guide/changelog' },
        ],
      },
      { text: 'API', link: '/api/' },
      { text: 'GitHub', link: 'https://github.com/szaghi/MOTIOn' },
    ],
    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',    link: '/guide/' },
            { text: 'Features', link: '/guide/features' },
          ],
        },
        {
          text: 'Getting Started',
          items: [
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Usage',        link: '/guide/usage' },
            { text: 'API Reference', link: '/guide/api-reference' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Contributing',      link: '/guide/contributing' },
            { text: 'Coverage Analysis', link: '/guide/coverage-analysis' },
            { text: 'Changelog',         link: '/guide/changelog' },
          ],
        },
      ],
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],
    },
    search: {
      provider: 'local',
    },
  },
  mermaid: {},
  vite: {
    build: {
      target: 'es2022',
    },
    optimizeDeps: {
      include: ['mermaid'],
      exclude: ['debug'],
    },
  },
})
