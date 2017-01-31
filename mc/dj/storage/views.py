import json
from django.http import JsonResponse

from . import utils as storage_utils

def post_data(request):
    data = request.POST.get('data', None)
    if isinstance(data, str): data = data.encode('utf-8')
    serialized_params = request.POST.get('params', None)
    params = _deserialize_params(serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    updated_params = backend.post_data(data=data, params=params)
    return JsonResponse({'params': updated_params})

def get_data(request):
    serialized_params = request.GET.get('params', None)
    params = _deserialize_params(serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    data = backend.get_data(params=params)
    if isinstance(data, bytes): data = data.decode('utf-8')
    return JsonResponse({'data': data})

def _deserialize_params(serialized_params=None):
    if serialized_params is None: return None
    return json.loads(serialized_params)
