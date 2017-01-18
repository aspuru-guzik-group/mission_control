from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets

from .models import Mol 
from .serializers import MolSerializer


class MolViewSet(viewsets.ModelViewSet):
    queryset = Mol.objects.all()
    serializer_class = MolSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_fields = ('uuid',)

@require_http_methods(["GET"])
def counts(request):
    counts = {
        'Mol': Mol.objects.count()
    }
    return JsonResponse(counts)
