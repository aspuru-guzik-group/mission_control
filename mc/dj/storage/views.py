from django.http import JsonResponse

from . import utils as storage_utils

def post_data(request):
    data = request.POST.get('data', None)
    params = request.POST.get('params', None)
    backend = storage_utils.get_storage_backend(params=params)
    storage_uri = backend.post(data=data, params=params)
    return JsonResponse({'storage_uri': storage_uri})

def get_data(request):
    params = request.GET.get('params', None)
    backend = storage_utils.get_storage_backend(params=params)
    return JsonResponse({'data': backend.get(params=params)})
