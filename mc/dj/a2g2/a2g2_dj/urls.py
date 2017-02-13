from django.conf.urls import url, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'chemthings', views.ChemThingViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^counts/', views.counts, name='counts'),
    url(r'^flush/', views.flush, name='flush'),
]
