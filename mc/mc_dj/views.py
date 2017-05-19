import json
from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from django_filters.rest_framework import DjangoFilterBackend, FilterSet
from rest_framework import viewsets
from rest_framework.decorators import detail_route

from mc_dj.constants import JobStatuses
from mc_dj import models as _models
from mc_dj import serializers as _serializers
from mc_dj.utils import queue_utils as _queue_utils


class FlowFilter(FilterSet):
    class Meta:
        model = _models.Flow
        fields = ['uuid', 'status']

    @property
    def qs(self):
        _qs = super().qs
        if self.request_has_tickable_param(self.request):
            _qs = self.filter_for_tickable(qs=_qs)
        return _qs

    def request_has_tickable_param(self, request=None):
        return hasattr(self.request, 'GET') and ('tickable' in request.GET)

    def filter_for_tickable(self, qs=None):
        tickable_statuses = [s.name
                             for s in _models.FlowStatuses.tickable_statuses]
        modified_qs = qs.filter(status__in=tickable_statuses)
        return modified_qs

class FlowViewSet(viewsets.ModelViewSet):
    queryset = _models.Flow.objects.all()
    serializer_class = _serializers.FlowSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_class = FlowFilter

@require_http_methods(["POST"])
@csrf_exempt
def claim_flows(request):
    result = {}
    post_data = json.loads(request.body.decode())
    uuids = post_data.get('uuids', [])
    if uuids:
        flows = _models.Flow.objects.filter(uuid__in=uuids)
        for flow in flows:
            if flow.claimed:
                result[flow.uuid] = None
            else:
                flow.claimed = True
                flow.save()
                result[flow.uuid] = _serializers.FlowSerializer(flow).data
    return JsonResponse(result)

class JobViewSet(viewsets.ModelViewSet):
    queryset = _models.Job.objects.all()
    serializer_class = _serializers.JobSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_fields = ('status', 'uuid',)

@require_http_methods(["POST"])
@csrf_exempt
def claim_jobs(request):
    result = {}
    post_data = json.loads(request.body.decode())
    uuids = post_data.get('uuids', [])
    if uuids:
        jobs = _models.Job.objects.filter(uuid__in=uuids)
        for job in jobs:
            if job.status == JobStatuses.PENDING.name:
                job.status = JobStatuses.RUNNING.name
                job.save()
                result[job.uuid] = _serializers.JobSerializer(job).data
            else:
                result[job.uuid] = None
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
    serializer_class = _serializers.QueueSerializer

    @detail_route(methods=['get', 'post'])
    def claim_items(self, request, pk=None):
        queue = self.get_object()
        models = {}
        serializers = {}
        for item_type in ['Job', 'Flow']:
            models[item_type] = getattr(_models, item_type)
            serializers[item_type] = getattr(_serializers,
                                             '%sSerializer' % item_type)
        claimed_items = _queue_utils.claim_queue_items(queue=queue,
                                                       models=models)
        serialized_items = _queue_utils.serialize_queue_items(
            queue=queue, items=claimed_items, serializers=serializers)
        result = {'items': serialized_items}
        return JsonResponse(result)
