import json
import urllib.parse

from django.http import HttpResponse, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.urls import reverse_lazy

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
    return JsonResponse({
        'params': updated_params,
        'url': get_url_for_params(params=updated_params)
    })

def get_url_for_params(params=None):
    url = "{path}?{query_str}".format(
        path=reverse_lazy('storage:storage-get_data'),
        query_str=urllib.parse.urlencode({
            'params': _serialize_params(params=params)
        })
    )
    return url

@csrf_exempt
def get_data(request):
    serialized_params = request.GET.get('params', None)
    params = _deserialize_params(serialized_params=serialized_params)
    backend = storage_utils.get_storage_backend(params=params)
    data = backend.get_data(params=params)
    return HttpResponse(data)

def _serialize_params(params=None):
    params or {}
    return json.dumps(params)

def _deserialize_params(serialized_params=None):
    if serialized_params is None: return None
    return json.loads(serialized_params)
