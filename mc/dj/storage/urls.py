from django.conf.urls import url
from . import views as _views

urlpatterns = (
    url(r'^post/$', _views.post_data),
    url(r'^get/$', _views.get_data),
)
