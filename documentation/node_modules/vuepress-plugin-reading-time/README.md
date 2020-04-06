# vuepress-plugin-reading-time

> reading time plugin to display how long a post takes to read

```bash
yarn add -D vuepress-plugin-reading-time
# OR npm install -D vuepress-plugin-reading-time
```

## Usage

```javascript
module.exports = {
  plugins: ['vuepress-plugin-reading-time']
}
```

Puts `reading-time` data into `$page` data so you can access like

```js
console.log($page.readingTime)
```

Example output
```js
{
  text: '1 min read',
  minutes: 0.08,
  time: 4800,
  words: 16
}
```

### Overide

You can override by specifying a readingTime object in frontmatter

```
---
title: My great post!
date: "2018-08-13T17:36:55.338Z"
readingTime: { text: "Not so long!", minutes: 3 }
---

Some content here...
```


## Options

### excludes

- Type: `Array<string>`
- Default: ``

Exclude pages by their path via a regular expression. This tests for both `path`
and `regularPath`.

Example:

```javascript
plugins: [
  ['vuepress-plugin-reading-time', {
    excludes: ['/about', '/tag/.*']
  }]
]
```
