from django.shortcuts import render
from django.http import HttpResponse

from .utils import task_registry

def create_task(request, *args, **kwargs):
    content = "\n".join([task._fw_name
                         for task in task_registry.tasks.values()])
    return HttpResponse(content)
