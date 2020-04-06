import Tabs from './components/Tabs.vue';
import TabPane from './components/Tab-pane.vue';

import './theme/tabs.scss';

export default function(ctx) {
  const { Vue } = ctx;
  Vue.component('Tabs', Tabs);
  Vue.component('Tab', TabPane);
}
