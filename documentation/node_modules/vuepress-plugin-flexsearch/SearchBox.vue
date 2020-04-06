<template>
  <div class="search-box">
    <input
      @input="query = $event.target.value"
      aria-label="Search"
      :value="query"
      :class="{ 'focused': focused }"
      :placeholder="placeholder"
      autocomplete="off"
      spellcheck="false"
      @focus="focused = true"
      @blur="focused = false"
      @keyup.enter="go(focusIndex)"
      @keyup.up="onUp"
      @keyup.down="onDown"
      ref="input"
    >
    <ul
      class="suggestions"
      v-if="showSuggestions"
      :class="{ 'align-right': alignRight }"
      @mouseleave="unfocus"
    >
      <li
        class="suggestion"
        v-for="(s, i) in suggestions"
        :key="i"
        :class="{ focused: i === focusIndex }"
        @mousedown="go(i)"
        @mouseenter="focus(i)"
      >
        <a
          :href="s.path"
          @click.prevent
        >
          <span
            class="page-title"
            v-html="s.title || s.path"
          ></span>
        </a>
      </li>
    </ul>
  </div>
</template>

<script>
import Flexsearch from "flexsearch";

/* global SEARCH_MAX_SUGGESTIONS, SEARCH_PATHS, SEARCH_HOTKEYS */
export default {
  data () {
    return {
      query: '',
      focused: false,
      focusIndex: 0,
      placeholder: undefined,
      index: null
    }
  },

  mounted () {
    this.placeholder = this.$site.themeConfig.searchPlaceholder || ''
    document.addEventListener('keydown', this.onHotkey)

    this.setupFlexSearch()
  },

  beforeDestroy () {
    document.removeEventListener('keydown', this.onHotkey)
  },

  computed: {
    showSuggestions () {
      return (
        this.focused
        && this.suggestions
        && this.suggestions.length
      )
    },

    suggestions () {
      const query = this.query.trim().toLowerCase()
      if (!query) {
        return
      }

      const result = this.index.search(query, 10).map(page => {
        return {
          ...page,
          title: this.getQuerySnippet(page)
        }
      })

      return result
    },

    // make suggestions align right when there are not enough items
    alignRight () {
      const navCount = (this.$site.themeConfig.nav || []).length
      const repo = this.$site.repo ? 1 : 0
      return navCount + repo <= 2
    }
  },

  methods: {
    getPageLocalePath (page) {
      for (const localePath in this.$site.locales || {}) {
        if (localePath !== '/' && page.path.indexOf(localePath) === 0) {
          return localePath
        }
      }
      return '/'
    },

    isSearchable (page) {
      let searchPaths = SEARCH_PATHS

      // all paths searchables
      if (searchPaths === null) { return true }

      searchPaths = Array.isArray(searchPaths) ? searchPaths : new Array(searchPaths)

      return searchPaths.filter(path => {
        return page.path.match(path)
      }).length > 0
    },

    onHotkey (event) {
      if (event.srcElement === document.body && SEARCH_HOTKEYS.includes(event.key)) {
        this.$refs.input.focus()
        event.preventDefault()
      }
    },

    onUp () {
      if (this.showSuggestions) {
        if (this.focusIndex > 0) {
          this.focusIndex--
        } else {
          this.focusIndex = this.suggestions.length - 1
        }
      }
    },

    onDown () {
      if (this.showSuggestions) {
        if (this.focusIndex < this.suggestions.length - 1) {
          this.focusIndex++
        } else {
          this.focusIndex = 0
        }
      }
    },

    go (i) {
      if (!this.showSuggestions) {
        return
      }
      const path = this.suggestions[i].path

      if (this.$route.path !== path) {
        this.$router.push(this.suggestions[i].path)
      };

      this.query = ''
      this.focusIndex = 0
    },

    focus (i) {
      this.focusIndex = i
    },

    unfocus () {
      this.focusIndex = -1
    },

    setupFlexSearch () {
      let defaultOptions = {
        encode: "extra",
        tokenize: "full",
        threshold: 1,
        resolution: 3
      }
      let options = this.$site.themeConfig.flexSearchOptions || defaultOptions

      options = {
        ...options,
        doc: {
          id: "key",
          field: ["title", "content"]
        }
      }
      this.index = new Flexsearch(options);
      const { pages } = this.$site;
      this.index.add(pages);
    },

    getQuerySnippet (page) {
      const queryPosition = page.content.toLowerCase().indexOf(this.query)
      const startIndex = queryPosition - 20 < 0 ? 0 : queryPosition - 20
      const endIndex = queryPosition + 30
      let querySnippet = page.content.slice(startIndex, endIndex)
        .toLowerCase()
        .replace(/[\W_]+/g, " ")

      const queryWords = this.query.split(' ').filter((v, i, a) => !!v && a.indexOf(v) === i); 
      for (const word of queryWords) {
        querySnippet = querySnippet.replace(new RegExp(word, 'gi'), `<strong class="text--primary">${word}</strong>`)
      }

      if (querySnippet) {
        return `<strong class="text--primary">${page.title}</strong> > .. ${querySnippet} ..`
      } else {
        return page.title
      }
    },

  }
}
</script>

<style lang="stylus">
.search-box {
  display: inline-block;
  position: relative;
  margin-right: 1rem;

  input {
    cursor: text;
    width: 10rem;
    height: 2rem;
    color: lighten($textColor, 25%);
    display: inline-block;
    border: 1px solid darken($borderColor, 10%);
    border-radius: 2rem;
    font-size: 0.9rem;
    line-height: 2rem;
    padding: 0 0.5rem 0 2rem;
    outline: none;
    transition: all 0.2s ease;
    background: #fff url('search.svg') 0.6rem 0.5rem no-repeat;
    background-size: 1rem;

    &:focus {
      cursor: auto;
      border-color: $accentColor;
    }
  }

  .suggestions {
    background: #fff;
    width: 20rem;
    position: absolute;
    top: 1.5rem;
    border: 1px solid darken($borderColor, 10%);
    border-radius: 6px;
    padding: 0.4rem;
    list-style-type: none;

    &.align-right {
      right: 0;
    }
  }

  .suggestion {
    line-height: 1.4;
    padding: 0.4rem 0.6rem;
    border-radius: 4px;
    cursor: pointer;

    a {
      white-space: normal;
      color: lighten($textColor, 35%);

      .page-title {
        font-weight: 600;
      }

      .header {
        font-size: 0.9em;
        margin-left: 0.25em;
      }
    }

    &.focused {
      background-color: #f3f4f5;

      a {
        color: $accentColor;
      }
    }
  }
}

@media (max-width: $MQNarrow) {
  .search-box {
    input {
      cursor: pointer;
      width: 0;
      border-color: transparent;
      position: relative;

      &:focus {
        cursor: text;
        left: 0;
        width: 10rem;
      }
    }
  }
}

// Match IE11
@media all and (-ms-high-contrast: none) {
  .search-box input {
    height: 2rem;
  }
}

@media (max-width: $MQNarrow) and (min-width: $MQMobile) {
  .search-box {
    .suggestions {
      left: 0;
    }
  }
}

@media (max-width: $MQMobile) {
  .search-box {
    margin-right: 0;

    input {
      left: 1rem;
    }

    .suggestions {
      right: 0;
    }
  }
}

@media (max-width: $MQMobileNarrow) {
  .search-box {
    .suggestions {
      width: calc(100vw - 4rem);
    }

    input:focus {
      width: 8rem;
    }
  }
}

.text--primary {
  color: $accentColor;
}
</style>
