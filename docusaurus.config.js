module.exports = {
  title: 'PopGen.jl',
  tagline: 'Population Genetics in Julia',
  url: 'https://pdimens.github.io',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  organizationName: 'pdimens', // Usually your GitHub org/user name.
  projectName: 'PopGen.jl', // Usually your repo name.
  themeConfig: {
    navbar: {
      prism: {
        defaultLanguage: 'julia',
        additionalLanguages: 'julia',
      },
      algolia: {
        apiKey: 'f6d1b3005e55708c6b33b80157908b05',
        indexName: 'popgen_jl',
        //appId: 'app-id', // Optional, if you run the DocSearch crawler on your own
        //algoliaOptions: {}, // Optional, if provided by Algolia
      },
      title: 'PopGen.jl',
      logo: {
        alt: 'PopGen.jl logo',
        src: 'img/logo_icon.png',
      },
      links: [
        {
          to: 'docs/',
          activeBasePath: 'docs',
          label: 'Docs',
          position: 'right',
        },
        {
          to: 'docs/getting_started/about',
          label: 'About',
          position: 'right',
        },        
        {to: 'blog', label: 'Blog', position: 'right'},
        {
          to: 'docs/getting_started/community',
          label: 'Get Involved',
          position: 'right',
        },
        {
          href: 'https://github.com/pdimens/popgen.jl',
          label: 'GitHub',
          position: 'right',
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
          homePageId: 'getting_started/install',
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/pdimens/popgen.jl/edit/documentation/',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/pdimens/popgen.jl/edit/edit/documentation/website/blog/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
