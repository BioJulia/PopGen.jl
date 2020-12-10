const remarkMath = require("remark-math");
const rehypeKatex = require("rehype-katex");
module.exports = {
  title: 'PopGen.jl',
  tagline: 'A package to perform a suite of population genetics analyses, written in Julia',
  url: 'https://pdimens.github.io',
  baseUrl: '/PopGen.jl/',
  favicon: 'img/favicon.ico',
  organizationName: 'pdimens', // Usually your GitHub org/user name.
  projectName: 'PopGen.jl', // Usually your repo name.
  themeConfig: {
    prism: {
      //defaultLanguage: 'julia',
      additionalLanguages: ['julia', 'r'],
      theme: require('prism-react-renderer/themes/github'),
      darkTheme: require('prism-react-renderer/themes/oceanicNext'),
    },
    algolia: {
      apiKey: 'f6d1b3005e55708c6b33b80157908b05',
      indexName: 'popgen_jl',
      //appId: 'app-id', // Optional, if you run the DocSearch crawler on your own
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
          to: 'docs/getting_started/about',
          label: 'About',
          position: 'left',
        },        
        {
          to: 'docs/getting_started/quickstart',
          label: 'Quickstart',
          position: 'right',
        },
        {to: 'blog', label: 'Blog', position: 'left'},
        {
        to: 'docs/latest',
        label: 'What\'s New',
        position: 'left',
        },
        {
          to: 'docs/getting_started/community',
          label: 'Get Involved',
          position: 'right',
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
          editUrl: 'https://github.com/pdimens/popgen.jl/edit/documentation/',
          showLastUpdateTime: true,
          remarkPlugins: [remarkMath],
          rehypePlugins: [[rehypeKatex, {strict: false}]],
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/pdimens/popgen.jl/edit/documentation/',
          feedOptions: {
            type: 'all',
            copyright: `Copyright © ${new Date().getFullYear()} PopGen.jl`,
          },
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
        sitemap: {
          cacheTime: 600 * 1000, // 600 sec - cache purge period
          changefreq: "weekly",
          priority: 0.5,
        },
      },
    ],
  ],
};
