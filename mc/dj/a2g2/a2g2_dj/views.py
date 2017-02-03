from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets

from .models import ChemThing 
from .serializers import ChemThingSerializer


class ChemThingViewSet(viewsets.ModelViewSet):
    queryset = ChemThing.objects.all()
    serializer_class = ChemThingSerializer
    filter_backends = (DjangoFilterBackend,)
    filter_fields = ('uuid',)

@require_http_methods(["GET"])
def counts(request):
    counts = {
        'ChemThing': ChemThing.objects.count()
    }
    return JsonResponse(counts)
