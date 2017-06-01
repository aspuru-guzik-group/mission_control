
from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
import django_filters as _dj_filters
import django_filters.rest_framework as _drf_filters
from rest_framework import viewsets
from rest_framework.decorators import list_route

from .models import ChemThing, a2g2_dj_models
from .serializers import ChemThingSerializer
from . import utils as _a2g2_dj_utils


class ChemThingFilter(_drf_filters.FilterSet):
    tag = _dj_filters.CharFilter(name="tags__name")
    class Meta:
        model = ChemThing
        fields = ['uuid', 'tag']

class ChemThingViewSet(viewsets.ModelViewSet):
    queryset = ChemThing.objects.all()
    serializer_class = ChemThingSerializer
    filter_backends = (_drf_filters.DjangoFilterBackend,)
    filter_class = ChemThingFilter

    @list_route(methods=['post'])
    def _bulk(self, request):
        result = _a2g2_dj_utils.process_serialized_chemthing_actions(
            serialized_chemthing_actions=request.body.decode())
        return JsonResponse(result, safe=False)

@require_http_methods(["GET"])
def counts(request):
    counts = {'ChemThing': ChemThing.objects.count()}
    return JsonResponse(counts)

@require_http_methods(["GET"])
def flush(request):
    flush_results = {}
    for model in a2g2_dj_models:
        model.objects.all().delete()
        flush_results[model.__name__] = 'flushed'
    return JsonResponse(flush_results)
