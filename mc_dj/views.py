import json
from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from django_filters.rest_framework import DjangoFilterBackend, FilterSet
from rest_framework import viewsets
from rest_framework.decorators import detail_route

from mc.mc_daos.dj_dao import DjDao
from .constants import Statuses
from . import models as _models
from . import serializers as _serializers

model_serializers = _serializers.get_model_serializers()

class FlowFilter(FilterSet):
    class Meta:
        model = _models.Flow
        fields = ['key', 'status']

    @property
    def qs(self):
        _qs = super().qs
        if self.request_has_tickable_param(self.request):
            _qs = self.filter_for_tickable(qs=_qs)
        return _qs

    def request_has_tickable_param(self, request=None):
        return hasattr(self.request, 'GET') and ('tickable' in request.GET)

    def filter_for_tickable(self, qs=None):
        tickable_statuses = [s for s in Statuses.tickable_statuses]
        modified_qs = qs.filter(status__in=tickable_statuses)
        return modified_qs

class FlowViewSet(viewsets.ModelViewSet):
    queryset = _models.Flow.objects.all()
    serializer_class = model_serializers['Flow']
    filter_backends = (DjangoFilterBackend,)
    filter_class = FlowFilter

@require_http_methods(["POST"])
@csrf_exempt
def claim_flows(request):
    result = {}
    post_data = json.loads(request.body.decode())
    keys = post_data.get('keys', [])
    if keys:
        flows = _models.Flow.objects.filter(key__in=keys)
        for flow in flows:
            if flow.claimed:
                result[flow.key] = None
            else:
                flow.claimed = True
                flow.save()
                result[flow.key] = _serializers.FlowSerializer(flow).data
    return JsonResponse(result)

class JobViewSet(viewsets.ModelViewSet):
    queryset = _models.Job.objects.all()
    serializer_class = model_serializers['Job']
    filter_backends = (DjangoFilterBackend,)
    filter_fields = ('status', 'key',)

@require_http_methods(["POST"])
@csrf_exempt
def claim_jobs(request):
    result = {}
    post_data = json.loads(request.body.decode())
    keys = post_data.get('keys', [])
    if keys:
        jobs = _models.Job.objects.filter(key__in=keys)
        for job in jobs:
            if job.status == Statuses.PENDING:
                job.status = Statuses.RUNNING
                job.save()
                result[job.key] = _serializers.JobSerializer(job).data
            else:
                result[job.key] = None
    return JsonResponse(result)

@require_http_methods(["GET"])
def flush(request):
    flush_results = {}
    for model in _models.mc_models:
        model.objects.all().delete()
        flush_results[model.__name__] = 'flushed'
    return JsonResponse(flush_results)

class QueueViewSet(viewsets.ModelViewSet):
    queryset = _models.Queue.objects.all()
    serializer_class = model_serializers['Queue']

    @detail_route(methods=['get', 'post'])
    def claim_items(self, request, pk=None):
        dao = get_dj_dao()
        result = {'items': dao.claim_queue_items(queue_key=pk)}
        return JsonResponse(result)

def get_dj_dao():
    item_types = ['Flow', 'Job', 'Queue']
    return DjDao(
        models={
            item_type: getattr(_models, item_type)
            for item_type in item_types
        },
        serializers={
            item_type: getattr(_serializers, item_type + 'Serializer')
            for item_type in item_types
        }
    )
