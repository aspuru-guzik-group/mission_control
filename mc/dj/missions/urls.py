from django.conf.urls import url, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'workflows', views.WorkflowViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^claim_workflows/', views.claim_workflows, name='claim_workflows'),
]
