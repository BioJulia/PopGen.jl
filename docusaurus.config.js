const remarkMath = require("remark-math");
const rehypeKatex = require("rehype-katex");
module.exports = {
  title: 'PopGen.jl',
  tagline: 'A package suite for population genetics analyses, written in Julia',
  url: 'https://BioJulia.dev',
  baseUrl: '/PopGen.jl/',
  favicon: 'img/favicon.ico',
  organizationName: 'BioJulia', // Usually your GitHub org/user name.
  projectName: 'PopGen.jl', // Usually your repo name.
  trailingSlash: false,
  plugins: [
    'plugin-image-zoom'
  ],
  stylesheets: [
    "https://fonts.googleapis.com/icon?family=Material+Icons",
  ],
  themeConfig: {
    docs : {
      sidebar: {
        hideable: true,
        autoCollapseCategories: true
    }},
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
    zoomSelector: '.markdown :not(em) > img',
    prism: {
      //defaultLanguage: 'julia',
      additionalLanguages: ['julia', 'r'],
      theme: require('prism-react-renderer/themes/github'),
      darkTheme: require('prism-react-renderer/themes/oceanicNext'),
    },
    algolia: {
      apiKey: '0fe48b3e529dd3af5c37979964b2ef41',
      indexName: 'popgen_jl',
      appId: 'JASY30KF23', // Optional, if you run the DocSearch crawler on your own
      //algoliaOptions: {}, // Optional, if provided by Algolia
    },
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
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          // It is recommended to set document id as docs home page (`docs/` path).
          //homePageId: 'docs/',
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl: 'https://github.com/BioJulia/PopGen.jl/edit/documentation/',
          showLastUpdateTime: true,
          remarkPlugins: [remarkMath],
          rehypePlugins: [[rehypeKatex, {strict: false}]],
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/BioJulia/PopGen.jl/edit/documentation/',
          feedOptions: {
            type: 'all',
            copyright: `Copyright © ${new Date().getFullYear()} PopGen.jl`,
          },
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
        sitemap: {
          //cacheTime: 600 * 1000, // 600 sec - cache purge period
          changefreq: "weekly",
          priority: 0.5,
        },
      },
    ],
  ],
};
