import json
from django.http import JsonResponse

from . import utils as storage_utils

def post_data(request):
    data = request.POST.get('data', None)
    serialized_params = request.POST.get('params', None)
    params = _deserialize_params(serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    storage_uri = backend.post_data(data=data, params=params)
    return JsonResponse({'storage_uri': storage_uri})

def get_data(request):
    serialized_params = request.GET.get('params', None)
    params = _deserialize_params(serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    return JsonResponse({'data': backend.get_data(params=params)})

def _deserialize_params(serialized_params=None):
    if serialized_params is None: return None
    return json.loads(serialized_params)
