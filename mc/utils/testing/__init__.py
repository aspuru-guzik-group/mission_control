from contextlib import contextmanager
import difflib
import filecmp
import io
import json
import os
import re
import sys
from unittest.mock import Mock


@contextmanager
def capture():
    with capture_std_streams(streams=['stdout']) as captured_streams:
        yield captured_streams['stdout']

def capture_std_streams(streams=None):
    stream_keys = streams or ['stdout']
    stream_cfgs = {}
    for stream_key in stream_keys:
        def _get_setter(stream_key_):
            return lambda new_stream: setattr(sys, stream_key_, new_stream)
        stream_cfgs[stream_key] = {'original': getattr(sys, stream_key),
                                   'setter': _get_setter(stream_key)}
    return capture_streams(stream_cfgs=stream_cfgs)

@contextmanager
def capture_streams(stream_cfgs=None):
    stream_cfgs = stream_cfgs or {}
    originals = {key: stream_cfg['original']
                 for key, stream_cfg in stream_cfgs.items()}
    try:
        tmp_streams = {}
        for key, stream_cfg in stream_cfgs.items():
            tmp_streams[key] = io.StringIO()
            stream_cfg['setter'](tmp_streams[key])
        yield tmp_streams
        for tmp_stream in tmp_streams.values(): tmp_stream.seek(0)
    finally:
        for key, original_stream in originals.items():
            stream_cfgs[key]['setter'](original_stream)

def assert_dirs_equal(test=None, left=None, right=None, ignore_patterns=None,
                      json_patterns=None):
    ignore_patterns = ignore_patterns or []
    json_patterns = json_patterns or [r'\.json$']
    cmp_result = filecmp.dircmp(left, right)
    test.assertEqual(set(cmp_result.left_only), set())
    test.assertEqual(set(cmp_result.right_only), set())
    for differing_file in cmp_result.diff_files:
        for ignore_pattern in ignore_patterns:
            if re.search(ignore_pattern, differing_file): continue
        left_file = os.path.join(left, differing_file)
        right_file = os.path.join(right, differing_file)
        check_as_json = False
        for json_pattern in json_patterns:
            if re.search(json_pattern, differing_file):
                check_as_json = True
                break
        if check_as_json:
            assert_json_files_equal(test=test, left_file=left_file,
                                    right_file=right_file)
        else:
            test.fail("Files differ:\n%s" % (
                diff(left_file=left_file, right_file=right_file)))

def assert_json_files_equal(test=None, left_file=None, right_file=None):
    test.assertEqual(json.load(open(left_file)), json.load(open(right_file)))

def diff(left_file=None, right_file=None):
    diff_lines = difflib.unified_diff(
        open(left_file).readlines(), open(right_file).readlines(),
        fromfile=left_file, tofile=right_file)
    return "".join([line for line in diff_lines])


def patch_request_client(request_client=None,
                         methods_to_patch=['get', 'patch', 'put', 'post']):
    for method_name in methods_to_patch:
        patch_client_method(client=request_client, method_name=method_name)
    request_client.raise_for_status = Mock()

def patch_client_method(client=None, method_name=None):
    orig_method = getattr(client, method_name)
    def patched_method(*args, data=None, files=None, **kwargs):
        args = list(args) or []
        if 'json' in kwargs:
            if method_name is not 'get':
                data = json.dumps(kwargs['json'])
                kwargs = {**kwargs, 'content_type': 'application/json'}
        else:
            data = data or {}
            if 'params' in kwargs: data.update(kwargs['params'])
            for name, requests_file in (files or []):
                data[name] = _requests_file_to_file_handle(requests_file)
        if data: args.append(data)
        response =  orig_method(*args, **kwargs)
        def raise_for_status():
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: ", response)
        response.raise_for_status = raise_for_status
        return response
    setattr(client, method_name, patched_method)

def _requests_file_to_file_handle(requests_file=None):
    name, file_data, content_type = requests_file
    if isinstance(file_data, bytes):
        file_handle = io.BytesIO(file_data)
    else:
        file_handle = io.StringIO(file_data)
    file_handle.content_type = content_type
    return file_handle
