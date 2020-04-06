'use strict';

var version = require('../version');
var getUA = function getUA() {
  var comments = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : [];

  var ua = 'LeanCloud-JS-SDK/' + version;
  if (comments.length) {
    ua += ' (' + comments.join('; ') + ')';
  }
  return ua;
};
module.exports = getUA;