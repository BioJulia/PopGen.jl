'use strict';

var _stringify = require('babel-runtime/core-js/json/stringify');

var _stringify2 = _interopRequireDefault(_stringify);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var _require = require('../adapter'),
    getAdapter = _require.getAdapter;

var debug = require('debug')('leancloud:qiniu');

module.exports = function (uploadInfo, data, file) {
  var saveOptions = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : {};

  file.attributes.url = uploadInfo.url;
  file._bucket = uploadInfo.bucket;
  file.id = uploadInfo.objectId;
  // Get the uptoken to upload files to qiniu.
  var uptoken = uploadInfo.token;
  var url = uploadInfo.upload_url || 'https://upload.qiniup.com';
  var fileFormData = {
    field: 'file',
    data: data,
    name: file.attributes.name
  };
  var options = {
    headers: file._uploadHeaders,
    data: {
      name: file.attributes.name,
      key: uploadInfo.key || file._qiniu_key,
      token: uptoken
    },
    onprogress: saveOptions.onprogress
  };
  debug('url: %s, file: %o, options: %o', url, fileFormData, options);
  var upload = getAdapter('upload');
  return upload(url, fileFormData, options).then(function (response) {
    debug(response.status, response.data);
    if (response.ok === false) {
      var message = response.status;
      if (response.data) {
        if (response.data.error) {
          message = response.data.error;
        } else {
          message = (0, _stringify2.default)(response.data);
        }
      }
      var error = new Error(message);
      error.response = response;
      throw error;
    }
    return file;
  }, function (error) {
    var response = error.response;

    if (response) {
      debug(response.status, response.data);
      error.statusCode = response.status;
      error.response = response.data;
    }
    throw error;
  });
};