module.exports = {
    base: '/PopGen.jl/',
    //theme: 'vuepress-theme-cool',
    title: 'PopGen.jl',
    description: 'Population Genetics in Julia',
    plugins: [
        'flexsearch',
        'vuepress-plugin-element-tabs',
        'vuepress-plugin-smooth-scroll',
        '@vuepress/plugin-back-to-top',
        '@vuepress/active-header-links',
        '@vuepress/plugin-nprogress',
        '@vuepress/medium-zoom',
        '@vuepress/last-updated',
        [
         'vuepress-plugin-mathjax', {
            target: 'svg',
            macros: {
                '*': '\\times',
            },
        },
      ]
    ],
    themeConfig: {
        logo: '/images/logo_icon.png',
        flexSearchOptions: {
          },
        nav: [
            { text: 'Home', link: '/' },
            { text: 'About', link: '/getting_started/about' },
            { text: 'Docs', link: '/getting_started/install' },
            { text: 'Contribute', link: '/community' },
            { text: 'GitHub', link: 'https://github.com/pdimens/PopGen.jl' },
        ],
        sidebar: [
            {
            title: 'Getting Started',   // required
            //collapsable: true, // optional, defaults to true
            //sidebarDepth: 1,    // optional, defaults to 1
            children: [
                '/getting_started/install',
                '/getting_started/julia_primer',
                '/getting_started/comparison',
                '/getting_started/popobj_type',
                '/getting_started/other_types',
            ]
        },
        {
            title: 'Importing Data',   // required
            //collapsable: true, // optional, defaults to true
            sidebarDepth: 2,    // optional, defaults to 1
            children: [
                '/guide/io/file_import',
                '/guide/io/delimited',
                '/guide/io/genepop',
                '/guide/io/variantcall',
                '/guide/io/datasets',
            ]
        },
        {
            title: 'Tutorials',
            //collapsable: true,
            //sidebarDepth: 1,
            children: [
                '/tutorials/manipulate',
                '/tutorials/accessing_popdata',
                '/tutorials/view_and_sort',
                '/tutorials/location_and_pop',
                '/tutorials/exclusion'
            ]
        },
        {
            title: 'Analyses',
            //collapsable: true,
            //sidebarDepth: 1,
            children: [
                '/analyses/hardyweinberg',
                '/analyses/relatedness',
            ]
        }
        ],
        docsRepo: 'pdimens/popgen.jl',
        // if your docs are not at the root of the repo:
        docsDir: 'docs',
        // if your docs are in a specific branch (defaults to 'master'):
        docsBranch: 'gh-pages',
        // defaults to false, set to true to enable
        editLinks: true,
        // custom text for edit link. Defaults to "Edit this page"
        editLinkText: 'Help us improve this page!',
        searchPlaceholder: 'Search the docs...'
    }
}