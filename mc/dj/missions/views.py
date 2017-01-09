from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django_filters.rest_framework import DjangoFilterBackend, FilterSet
from rest_framework import viewsets

from .models import Flow, FlowStatuses
from .serializers import FlowSerializer

class FlowFilter(FilterSet):
    class Meta:
        model = Flow
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
        tickable_statuses = [s.name for s in FlowStatuses.tickable_statuses]
        modified_qs = qs.filter(status__in=tickable_statuses)
        return modified_qs

class FlowViewSet(viewsets.ModelViewSet):
    queryset = Flow.objects.all()
    serializer_class = FlowSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_class = FlowFilter

@require_http_methods(["POST"])
def claim_flows(request):
    result = {}
    csv_uuids = request.POST['uuids']
    if csv_uuids:
        uuids = csv_uuids.split(',')
        flows = Flow.objects.filter(uuid__in=uuids)
        for flow in flows:
            if flow.claimed:
                result[flow.uuid] = None
            else:
                flow.claimed = True
                flow.save()
                result[flow.uuid] = FlowSerializer(flow).data
    return JsonResponse(result)
