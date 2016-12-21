from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django_filters.rest_framework import DjangoFilterBackend, FilterSet
from rest_framework import viewsets

from .models import Workflow, WorkflowStatuses
from .serializers import WorkflowSerializer

class WorkflowFilter(FilterSet):
    class Meta:
        model = Workflow
        fields = ['status']

    @property
    def qs(self):
        _qs = super().qs
        if self.request_has_tickable_param(self.request):
            _qs = self.filter_for_tickable(qs=_qs)
        return _qs

    def request_has_tickable_param(self, request=None):
        return hasattr(self.request, 'GET') and ('tickable' in request.GET)

    def filter_for_tickable(self, qs=None):
        tickable_statuses = [s.name for s in WorkflowStatuses.tickable_statuses]
        modified_qs = qs.filter(status__in=tickable_statuses)
        return modified_qs

class WorkflowViewSet(viewsets.ModelViewSet):
    queryset = Workflow.objects.all()
    serializer_class = WorkflowSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_class = WorkflowFilter

@require_http_methods(["POST"])
def claim_workflows(request):
    result = {}
    csv_uuids = request.POST['uuids']
    if csv_uuids:
        uuids = csv_uuids.split(',')
        workflows = Workflow.objects.filter(uuid__in=uuids)
        for workflow in workflows:
            if workflow.claimed:
                result[workflow.uuid] = None
            else:
                workflow.claimed = True
                workflow.save()
                result[workflow.uuid] = WorkflowSerializer(workflow).data
    return JsonResponse(result)
