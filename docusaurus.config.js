module.exports = {
  title: 'PopGen.jl',
  tagline: 'Population Genetics in Julia',
  url: 'https://your-docusaurus-test-site.com',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  organizationName: 'pdimens', // Usually your GitHub org/user name.
  projectName: 'PopGen.jl', // Usually your repo name.
  themeConfig: {
    navbar: {
      prism: {
        additionalLanguages: ['julia'],
      },
      title: 'PopGen.jl',
      logo: {
        alt: 'PopGen.jl Logo',
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
            'https://github.com/facebook/docusaurus/edit/master/website/',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/facebook/docusaurus/edit/master/website/blog/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
