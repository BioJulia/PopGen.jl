import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

const config: Config = {
  title: 'PopGen.jl',
  tagline: 'A package suite for population genetics analyses, written in Julia',
  url: 'https://BioJulia.dev',
  baseUrl: '/PopGen.jl/',
  favicon: 'img/favicon.ico',
  organizationName: 'BioJulia', // Usually your GitHub org/user name.
  projectName: 'PopGen.jl', // Usually your repo name.
  trailingSlash: false,

  // Required as of v3: what to do with broken links (was implicit "warn" before)
  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  plugins: [
    [
      '@cmfcmf/docusaurus-search-local',
      {
        indexDocSidebarParentCategories: 0,
      },
    ],
    'docusaurus-plugin-image-zoom',
  ],
  stylesheets: [
    'https://fonts.googleapis.com/icon?family=Material+Icons',
  ],
  themeConfig: {
    docs: {
      sidebar: {
        hideable: true,
        autoCollapseCategories: true,
      },
    },
    colorMode: {
      defaultMode: 'light',
      disableSwitch: false,
      respectPrefersColorScheme: true,
    },
    announcementBar: {
      id: 'supportus',
      content:
        '🔵 Like PopGen.jl? Give it a ⭐ on <a target="_blank" rel="noopener noreferrer" href="https://github.com/BioJulia/PopGen.jl">GitHub!</a>! 🟣',
    },
    // renamed themeConfig key for docusaurus-plugin-image-zoom v3
    zoom: {
      selector: '.markdown :not(em) > img',
      background: {
        light: 'rgb(255, 255, 255)',
        dark: 'rgb(50, 50, 50)',
      },
    },
    prism: {
      additionalLanguages: ['julia', 'r'],
      theme: prismThemes.github,
      darkTheme: prismThemes.oceanicNext,
    },
    // algolia: {
    //  apiKey: '0fe48b3e529dd3af5c37979964b2ef41',
    //  indexName: 'popgen_jl',
    //  appId: 'JASY30KF23', // Optional, if you run the DocSearch crawler on your own
    //  algoliaOptions: {}, // Optional, if provided by Algolia
    //  },
    navbar: {
      hideOnScroll: true,
      title: 'PopGen.jl',
      logo: {
        alt: 'PopGen.jl logo',
        src: 'img/logo_icon.png',
      },
      items: [
        {
          to: 'docs/',
          activeBasePath: 'docs',
          label: 'Docs',
          position: 'left',
        },
        {
          to: 'docs/gettingstarted/quickstart',
          label: 'Quickstart',
          position: 'left',
        },
        {
          to: 'docs/gettingstarted/about',
          label: 'About',
          position: 'left',
        },
        {to: 'blog', label: 'Blog', position: 'left'},
        {
          to: 'docs/latest',
          label: 'News',
          position: 'left',
        },
        {
          to: 'docs/gettingstarted/community',
          label: 'Get Involved',
          position: 'left',
        },
        {
          href: 'https://github.com/pdimens/popgen.jl',
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
      ],
    },
    footer: {
      style: 'dark',
      copyright: `Copyright © ${new Date().getFullYear()} Pavel Dimens & Jason Selwyn. Built with Docusaurus. Rawr!`,
    },
  } satisfies Preset.ThemeConfig,
  presets: [
    [
      'classic',
      {
        docs: {
          sidebarPath: './sidebars.js',
          editUrl: 'https://github.com/BioJulia/PopGen.jl/edit/documentation/',
          showLastUpdateTime: true,
          remarkPlugins: [remarkMath],
          rehypePlugins: [[rehypeKatex, {strict: false}]],
        },
        blog: {
          showReadingTime: true,
          editUrl: 'https://github.com/BioJulia/PopGen.jl/edit/documentation/',
          feedOptions: {
            type: 'all',
            limit: false, // keep old "unlimited feed" behavior (v3 default is 20)
            copyright: `Copyright © ${new Date().getFullYear()} PopGen.jl`,
          },
        },
        theme: {
          customCss: './src/css/custom.css',
        },
        sitemap: {
          changefreq: 'weekly',
          priority: 0.5,
        },
      } satisfies Preset.Options,
    ],
  ],
};

export default config;