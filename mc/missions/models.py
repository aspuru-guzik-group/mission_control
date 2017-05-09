import enum
import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel

from jsonfield import JSONField

from .constants import JobStatuses


def str_uuid4(): return str(uuid.uuid4())

def uuid_model_str(instance):
    return '<{class_name}: {uuid}>'.format(
        class_name=instance.__class__.__name__,
        uuid=instance.uuid)

missions_models = []

class Mission(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    label = models.CharField(null=True, max_length=1024)

    def __str__(self):
        return uuid_model_str(self)

missions_models.append(Mission)

class FlowStatuses(enum.Enum):
    PENDING = {'label': 'pending'}
    RUNNING = {'label': 'running'}
    COMPLETED = {'label': 'completed'}
    FAILED = {'label': 'failed'}
FlowStatuses.tickable_statuses = [getattr(FlowStatuses, status)
                                      for status in ['PENDING', 'RUNNING']]

class Flow(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    label = models.CharField(max_length=256, blank=True, null=True)
    serialization = models.TextField(null=True)
    mission = models.ForeignKey('Mission', null=True, on_delete=models.CASCADE)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in FlowStatuses],
                              default=FlowStatuses.PENDING.name)
    claimed = models.NullBooleanField(null=True, default=False)

    def __str__(self): return uuid_model_str(self)

missions_models.append(Flow)

class Job(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    label = models.CharField(max_length=256, blank=True, null=True)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in JobStatuses],
                              default=JobStatuses.PENDING.name)
    claimed = models.NullBooleanField(null=True, default=False)
    job_spec = JSONField(default=dict, null=True, blank=True)
    data = models.TextField(null=True, default='{}')

missions_models.append(Job)

class Queue(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    label = models.CharField(max_length=256, blank=True, null=True)
    queue_spec = JSONField(default=dict, blank=True, null=True)

missions_models.append(Queue)
