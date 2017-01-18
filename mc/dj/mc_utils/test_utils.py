import json
from unittest.mock import Mock


def patch_request_client(request_client=None, json_methods=['patch', 'post']):
    for method_name in json_methods:
        patch_client_method_to_use_json(client=request_client,
                                        method_name=method_name)
    request_client.raise_for_status = Mock()

def patch_client_method_to_use_json(client=None, method_name=None):
    orig_method = getattr(client, method_name)
    def patched_method(path, data=None, **kwargs):
        if method_name == 'get':
            response = orig_method(path, content_type='application/json',
                                   **kwargs)
        else:
            response =  orig_method(path, json.dumps(data), 
                                    content_type='application/json', **kwargs)
        def raise_for_status():
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: ", response)
        response.raise_for_status = raise_for_status
        return response
    setattr(client, method_name, patched_method)

