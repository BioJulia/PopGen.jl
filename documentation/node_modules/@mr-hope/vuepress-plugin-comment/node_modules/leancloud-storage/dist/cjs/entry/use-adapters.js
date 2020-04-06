'use strict';

var adapters = require('@leancloud/platform-adapters-node');
var getUA = require('../ua');
var comments = (process.env.CLIENT_PLATFORM ? [process.env.CLIENT_PLATFORM] : []).concat(require('../ua/comments'));

module.exports = function (AV) {
  AV.setAdapters(adapters);
  AV._sharedConfig.userAgent = getUA(comments);
  return AV;
};