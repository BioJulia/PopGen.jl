'use strict';

var _promise = require('babel-runtime/core-js/promise');

var _promise2 = _interopRequireDefault(_promise);

var _getPrototypeOf = require('babel-runtime/core-js/object/get-prototype-of');

var _getPrototypeOf2 = _interopRequireDefault(_getPrototypeOf);

var _stringify = require('babel-runtime/core-js/json/stringify');

var _stringify2 = _interopRequireDefault(_stringify);

var _typeof2 = require('babel-runtime/helpers/typeof');

var _typeof3 = _interopRequireDefault(_typeof2);

var _keys = require('babel-runtime/core-js/object/keys');

var _keys2 = _interopRequireDefault(_keys);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var _ = require('underscore');

var _require = require('promise-timeout'),
    timeout = _require.timeout;

var debug = require('debug');
var debugRequest = debug('leancloud:request');
var debugRequestError = debug('leancloud:request:error');

var _require2 = require('../adapter'),
    getAdapter = _require2.getAdapter;

var requestsCount = 0;

var ajax = function ajax(_ref) {
  var method = _ref.method,
      url = _ref.url,
      query = _ref.query,
      data = _ref.data,
      _ref$headers = _ref.headers,
      headers = _ref$headers === undefined ? {} : _ref$headers,
      time = _ref.timeout,
      onprogress = _ref.onprogress;

  if (query) {
    var queryString = (0, _keys2.default)(query).map(function (key) {
      var value = query[key];
      if (value === undefined) return undefined;
      var v = (typeof value === 'undefined' ? 'undefined' : (0, _typeof3.default)(value)) === 'object' ? (0, _stringify2.default)(value) : value;
      return encodeURIComponent(key) + '=' + encodeURIComponent(v);
    }).filter(function (qs) {
      return qs;
    }).join('&');
    url = url + '?' + queryString;
  }

  var count = requestsCount++;
  debugRequest('request(%d) %s %s %o %o %o', count, method, url, query, data, headers);

  var request = getAdapter('request');
  var promise = request(url, { method: method, headers: headers, data: data, onprogress: onprogress }).then(function (response) {
    debugRequest('response(%d) %d %O %o', count, response.status, response.data || response.text, response.header);
    if (response.ok === false) {
      var error = new Error();
      error.response = response;
      throw error;
    }
    return response.data;
  }).catch(function (error) {
    if (error.response) {
      if (!debug.enabled('leancloud:request')) {
        debugRequestError('request(%d) %s %s %o %o %o', count, method, url, query, data, headers);
      }
      debugRequestError('response(%d) %d %O %o', count, error.response.status, error.response.data || error.response.text, error.response.header);
      error.statusCode = error.response.status;
      error.responseText = error.response.text;
      error.response = error.response.data;
    }
    throw error;
  });
  return time ? timeout(promise, time) : promise;
};

// Helper function to check null or undefined.
var isNullOrUndefined = function isNullOrUndefined(x) {
  return _.isNull(x) || _.isUndefined(x);
};

var ensureArray = function ensureArray(target) {
  if (_.isArray(target)) {
    return target;
  }
  if (target === undefined || target === null) {
    return [];
  }
  return [target];
};

var transformFetchOptions = function transformFetchOptions() {
  var _ref2 = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {},
      keys = _ref2.keys,
      include = _ref2.include,
      includeACL = _ref2.includeACL;

  var fetchOptions = {};
  if (keys) {
    fetchOptions.keys = ensureArray(keys).join(',');
  }
  if (include) {
    fetchOptions.include = ensureArray(include).join(',');
  }
  if (includeACL) {
    fetchOptions.returnACL = includeACL;
  }
  return fetchOptions;
};

var getSessionToken = function getSessionToken(authOptions) {
  if (authOptions.sessionToken) {
    return authOptions.sessionToken;
  }
  if (authOptions.user && typeof authOptions.user.getSessionToken === 'function') {
    return authOptions.user.getSessionToken();
  }
};

var tap = function tap(interceptor) {
  return function (value) {
    return interceptor(value), value;
  };
};

// Shared empty constructor function to aid in prototype-chain creation.
var EmptyConstructor = function EmptyConstructor() {};

// Helper function to correctly set up the prototype chain, for subclasses.
// Similar to `goog.inherits`, but uses a hash of prototype properties and
// class properties to be extended.
var inherits = function inherits(parent, protoProps, staticProps) {
  var child;

  // The constructor function for the new subclass is either defined by you
  // (the "constructor" property in your `extend` definition), or defaulted
  // by us to simply call the parent's constructor.
  if (protoProps && protoProps.hasOwnProperty('constructor')) {
    child = protoProps.constructor;
  } else {
    /** @ignore */
    child = function child() {
      parent.apply(this, arguments);
    };
  }

  // Inherit class (static) properties from parent.
  _.extend(child, parent);

  // Set the prototype chain to inherit from `parent`, without calling
  // `parent`'s constructor function.
  EmptyConstructor.prototype = parent.prototype;
  child.prototype = new EmptyConstructor();

  // Add prototype properties (instance properties) to the subclass,
  // if supplied.
  if (protoProps) {
    _.extend(child.prototype, protoProps);
  }

  // Add static properties to the constructor function, if supplied.
  if (staticProps) {
    _.extend(child, staticProps);
  }

  // Correctly set child's `prototype.constructor`.
  child.prototype.constructor = child;

  // Set a convenience property in case the parent's prototype is
  // needed later.
  child.__super__ = parent.prototype;

  return child;
};

var parseDate = function parseDate(iso8601) {
  return new Date(iso8601);
};

var setValue = function setValue(target, key, value) {
  // '.' is not allowed in Class keys, escaping is not in concern now.
  var segs = key.split('.');
  var lastSeg = segs.pop();
  var currentTarget = target;
  segs.forEach(function (seg) {
    if (currentTarget[seg] === undefined) currentTarget[seg] = {};
    currentTarget = currentTarget[seg];
  });
  currentTarget[lastSeg] = value;
  return target;
};

var findValue = function findValue(target, key) {
  var segs = key.split('.');
  var firstSeg = segs[0];
  var lastSeg = segs.pop();
  var currentTarget = target;
  for (var i = 0; i < segs.length; i++) {
    currentTarget = currentTarget[segs[i]];
    if (currentTarget === undefined) {
      return [undefined, undefined, lastSeg];
    }
  }
  var value = currentTarget[lastSeg];
  return [value, currentTarget, lastSeg, firstSeg];
};

var isPlainObject = function isPlainObject(obj) {
  return _.isObject(obj) && (0, _getPrototypeOf2.default)(obj) === Object.prototype;
};

var continueWhile = function continueWhile(predicate, asyncFunction) {
  if (predicate()) {
    return asyncFunction().then(function () {
      return continueWhile(predicate, asyncFunction);
    });
  }
  return _promise2.default.resolve();
};

module.exports = {
  ajax: ajax,
  isNullOrUndefined: isNullOrUndefined,
  ensureArray: ensureArray,
  transformFetchOptions: transformFetchOptions,
  getSessionToken: getSessionToken,
  tap: tap,
  inherits: inherits,
  parseDate: parseDate,
  setValue: setValue,
  findValue: findValue,
  isPlainObject: isPlainObject,
  continueWhile: continueWhile
};