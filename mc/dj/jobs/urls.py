from django.conf.urls import url, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'jobs', views.JobViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^claim_jobs/', views.claim_jobs, name='claim_jobs'),
]
