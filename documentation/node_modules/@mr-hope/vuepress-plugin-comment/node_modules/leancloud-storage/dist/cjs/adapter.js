'use strict';

var _ = require('underscore');
var adapters = {};

var getAdapter = function getAdapter(name) {
  var adapter = adapters[name];
  if (adapter === undefined) {
    throw new Error(name + ' adapter is not configured');
  }
  return adapter;
};
var setAdapters = function setAdapters(newAdapters) {
  _.extend(adapters, newAdapters);
};

module.exports = {
  getAdapter: getAdapter,
  setAdapters: setAdapters
};