import json
from unittest.mock import Mock


def patch_request_client(request_client=None, json_methods=['patch', 'post']):
    for method_name in json_methods:
        patch_client_method_to_use_json(client=request_client,
                                        method_name=method_name)
    request_client.raise_for_status = Mock()

def patch_client_method_to_use_json(client=None, method_name=None):
    orig_method = getattr(client, method_name)
    def patched_method(*args, data=None, **kwargs):
        args = list(args) or []
        if kwargs.get('content_type', None) == 'application/json':
            if method_name is not 'get' and data: data = json.dumps(data)
        if data: args.append(data)
        response =  orig_method(*args, **kwargs)
        def raise_for_status():
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: ", response)
        response.raise_for_status = raise_for_status
        return response
    setattr(client, method_name, patched_method)

