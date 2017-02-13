import json
from django.http import HttpResponse, JsonResponse
from django.views.decorators.csrf import csrf_exempt

from . import utils as storage_utils

@csrf_exempt
def post_data(request):
    params = {}
    if 'params' in request.FILES:
        serialized_params = request.FILES['params'].read()
        if serialized_params:
            if isinstance(serialized_params, bytes):
                serialized_params = serialized_params.decode()
            params = json.loads(serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    updated_params = backend.post_data(data=request.FILES['data'],
                                       params=params)
    return JsonResponse({'params': updated_params})

@csrf_exempt
def get_data(request):
    serialized_params = request.GET.get('params', None)
    params = _deserialize_params(serialized_params=serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    data = backend.get_data(params=params)
    return HttpResponse(data)

def _deserialize_params(serialized_params=None):
    if serialized_params is None: return None
    return json.loads(serialized_params)
