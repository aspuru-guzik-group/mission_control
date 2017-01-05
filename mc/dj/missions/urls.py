from django.conf.urls import url, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'flows', views.FlowViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^claim_flows/', views.claim_flows, name='claim_flows'),
]
