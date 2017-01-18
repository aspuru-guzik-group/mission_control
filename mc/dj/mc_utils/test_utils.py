import json


def patch_request_client_to_use_json(client=None,
                                     methods_to_patch=['patch', 'post']):
    for method_name in methods_to_patch:
        patch_client_method_to_use_json(client=client, method_name=method_name)

def patch_client_method_to_use_json(client=None, method_name=None):
    orig_method = getattr(client, method_name)
    def patched_method(path, data=None, **kwargs):
        return orig_method(path, json.dumps(data), 
                           content_type='application/json', **kwargs)
    setattr(client, method_name, patched_method)
