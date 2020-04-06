module.exports = {
    theme: 'default-prefers-color-scheme',
    title: 'PopGen.jl',
    description: 'Population Genetics in Julia',
    plugins: [
        'flexsearch',
        'vuepress-plugin-element-tabs',
        'vuepress-plugin-smooth-scroll',
        '@vuepress/plugin-back-to-top',
        '@vuepress/active-header-links',
    //'vuepress-plugin-reading-time',
        '@vuepress/plugin-nprogress',
        ['@vuepress/medium-zoom', {
          selector: '.theme-default-content :not(a) > img'
        }],
        [
            'vuepress-plugin-clean-urls', {        
                normalSuffix: '',
                notFoundPath: '/404.html',},
    ],
    [
        'vuepress-plugin-mathjax',
        {
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
            // to override the default options you can see available options on https://github.com/nextapps-de/flexsearch
          },
        nav: [
            { text: 'Home', link: '/' },
            { text: 'About', link: '/guide/about' },
            { text: 'Docs', link: '/guide/' },
            { text: 'Contribute', link: '/community' },
            { text: 'GitHub', link: 'https://github.com/pdimens/PopGen.jl' },
        ],
        sidebar: [
            {
            title: 'Getting Started',   // required
            //path: '/guide/',      // optional, which should be a absolute path.
            collapsable: true, // optional, defaults to true
            sidebarDepth: 1,    // optional, defaults to 1
            children: [
                '/guide/',
                '/guide/install',
                '/guide/comparison',
                '/guide/popobj_type',
                '/guide/other_types',
            ]
        },
        {
            title: 'Importing Data',   // required
            //path: '/guide/',      // optional, which should be a absolute path.
            collapsable: true, // optional, defaults to true
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
            collapsable: true,
            sidebarDepth: 1,
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
            collapsable: true,
            sidebarDepth: 1,
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
        docsBranch: 'master',
        // defaults to false, set to true to enable
        editLinks: true,
        // custom text for edit link. Defaults to "Edit this page"
        editLinkText: 'Help us improve this page!',
        searchPlaceholder: 'Search the docs...'
    }
}