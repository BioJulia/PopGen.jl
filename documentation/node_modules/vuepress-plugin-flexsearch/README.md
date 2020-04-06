# vuepress-plugin-flexsearch

Next-Generation full text search library for Vuepress

> A good replacement of the default Vuepress search plugin.

## Installation

```bash
yarn add -D vuepress-plugin-flexsearch
// or npm install vuepress-plugin-flexsearch -D

```

## Usage

```js
// .vuepress/config.js

module.exports = {
    themeConfig: {
      flexSearchOptions: {
        // to override the default options you can see available options on https://github.com/nextapps-de/flexsearch
      }
    },
    plugins: [
      ['flexsearch'],
      // other plugins
    ]
}
```

## Credit

Thanks to [nextapps-de/flexsearch](https://github.com/nextapps-de/flexsearch)

## Contributions

PRs are welcome :heart:
