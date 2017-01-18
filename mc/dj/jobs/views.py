import json

from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets

from .models import Job, JobStatuses
from .serializers import JobSerializer


class JobViewSet(viewsets.ModelViewSet):
    queryset = Job.objects.all()
    serializer_class = JobSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_fields = ('status',)

@require_http_methods(["POST"])
def claim_jobs(request):
    result = {}
    post_data = json.loads(request.body.decode())
    uuids = post_data.get('uuids', [])
    if uuids:
        jobs = Job.objects.filter(uuid__in=uuids)
        for job in jobs:
            if job.status == JobStatuses.PENDING.name:
                job.status = JobStatuses.CLAIMED.name
                job.save()
                result[job.uuid] = JobSerializer(job).data
            else:
                result[job.uuid] = None
    return JsonResponse(result)
