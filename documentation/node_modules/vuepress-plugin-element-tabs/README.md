# [vuepress-plugin-element-tabs](https://superbiger.github.io/vuepress-plugin-tabs/)

<a href="https://www.npmjs.com/package/vuepress-plugin-element-tabs"><img src="https://img.shields.io/npm/v/vuepress-plugin-element-tabs.svg" alt="Version"></a>
<a href="https://www.npmjs.com/package/vuepress-plugin-element-tabs"><img src="https://img.shields.io/npm/l/vuepress-plugin-element-tabs.svg" alt="License"></a>
<img src="https://img.shields.io/badge/thanks-element-brightgreen.svg"/>

Vuepress plugin - markdown custom container to display content in tabs from [Element UI](https://github.com/ElemeFE/element)

## Docs
> https://superbiger.github.io/vuepress-plugin-tabs

## Install
> This plugin requires VuePress >= 1.0.0, for now you can try it via yarn add vuepress@next -D 

```shell
yarn add vuepress-plugin-element-tabs -D
```

```javascript
// .vuepress/config.js
module.exports = {
  plugins: [
    'vuepress-plugin-element-tabs'
  ]
}
```

## Usage

~~~ md
:::: tabs

::: tab title
__markdown content__
:::


::: tab javascript
``` javascript
() => {
  console.log('Javascript code example')
}
```
:::

::::

~~~

## Documents
> Accepted Value Like That
~~~md
:::: tabs type:board-card
::: tab title lazy
__markdown content__
:::
::::
~~~

### Tabs Attributes
|Attribute|Description|Type|Accepted Values|Default|
|:--|:--|:--|:--|:--|
|type|type of Tab|String|card/border-card|border-card|
|tab-position|position of tabs|String|top/right/bottom/left|top|
|stretch|whether width of tab automatically fits its container|Boolean|-|false|


### Tab Attributes
|Attribute|Description|Type|Accepted Values|Default|
|:--|:--|:--|:--|:--|
|label|title of the tab|String|-|-|
|lazy|whether Tab is lazily rendered|Boolean|-|false|

### Q&A
* How to get mouse position with canvas ?
```javascript

var canvas = document.getElementById('screen');
var mouse = getMouse(canvas);

function addEvent(obj, type, handle) {
  try {
    obj.addEventListener(type, handle, false);
  } catch (e) {
    try {
      obj.attachEvent("on" + type, handle);
    } catch (e) {
      obj["on" + type] = handle;
    }
  }
}

function getMouse(element) {
  var mouse = { x: 0, y: 0 };

  addEvent(element, "mousemove", function(e) {
    var x, y;
    var e = e || window.event;
    const box = element.getBoundingClientRect();
    x = e.clientX - box.x;
    y = e.clientY - box.y;
    mouse.x = x;
    mouse.y = y;
  });

  return mouse;
}
```

