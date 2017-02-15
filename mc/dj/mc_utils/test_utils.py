import io
import json
from unittest.mock import Mock


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
