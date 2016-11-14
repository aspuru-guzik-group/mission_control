from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^create/', views.create_task, name='task_create'),
]
