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
          position: 'left',
        },
        {to: 'blog', label: 'Blog', position: 'left'},
        {
          href: 'https://github.com/pdimens/popgen.jl',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Style Guide',
              to: 'docs/',
            },
            {
              label: 'Second Doc',
              to: 'docs/doc2/',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'PopGen.jl Slack',
              href: 'https://join.slack.com/t/popgenjl/shared_invite/zt-deam65n8-DuBs2z1oDtsbBuRplJW~Pg',
            },
            {
              label: 'Pavel\'s Twitter',
              href: 'https://twitter.com/PVDimens',
            },
            {
              label: 'Jason\'s Twitter',
              href: 'https://twitter.com/JasonSelwyn',
            }
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'Blog',
              to: 'blog',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/pdimens/popgen.jl',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Pavel Dimens & Jason Selwyn. Built with Docusaurus. Rawr!`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          // It is recommended to set document id as docs home page (`docs/` path).
          homePageId: 'doc1',
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
