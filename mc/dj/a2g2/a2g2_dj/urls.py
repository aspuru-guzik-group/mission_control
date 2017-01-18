from django.conf.urls import url, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'mols', views.MolViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^counts/', views.counts, name='counts'),
]
