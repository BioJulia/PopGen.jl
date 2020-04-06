'use strict';

var _promise = require('babel-runtime/core-js/promise');

var _promise2 = _interopRequireDefault(_promise);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var _require = require('./adapter'),
    getAdapter = _require.getAdapter;

var syncApiNames = ['getItem', 'setItem', 'removeItem', 'clear'];

var localStorage = {
  get async() {
    return getAdapter('storage').async;
  }
};

// wrap sync apis with async ones.
syncApiNames.forEach(function (apiName) {
  localStorage[apiName + 'Async'] = function () {
    var storage = getAdapter('storage');
    return _promise2.default.resolve(storage[apiName].apply(storage, arguments));
  };

  localStorage[apiName] = function () {
    var storage = getAdapter('storage');
    if (!storage.async) {
      return storage[apiName].apply(storage, arguments);
    }
    var error = new Error('Synchronous API [' + apiName + '] is not available in this runtime.');
    error.code = 'SYNC_API_NOT_AVAILABLE';
    throw error;
  };
});

module.exports = localStorage;