from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from rest_framework import viewsets

from .models import Job
from .serializers import JobSerializer


class JobViewSet(viewsets.ModelViewSet):
    queryset = Job.objects.all()
    serializer_class = JobSerializer

@require_http_methods(["POST"])
def claim_jobs(request):
    result = {}
    csv_uuids = request.POST['uuids']
    if csv_uuids:
        uuids = csv_uuids.split(',')
        jobs = Job.objects.filter(uuid__in=uuids)
        for job in jobs:
            if job.status is Job.STATUSES.PENDING.name:
                job.status = Job.STATUSES.CLAIMED.name
                job.save()
                result[job.uuid] = True
            else:
                result[job.uuid] = False
    return JsonResponse(result)
