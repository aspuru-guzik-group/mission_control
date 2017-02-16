import enum
import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import JSONField

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
    name = models.CharField(null=True, max_length=1024)

    def __str__(self):
        return uuid_model_str(self)

missions_models.append(Mission)

class FlowStatuses(enum.Enum):
    PENDING = {'label': 'pending'}
    RUNNING = {'label': 'running'}
    COMPLETED = {'label': 'completed'}
FlowStatuses.tickable_statuses = [getattr(FlowStatuses, status)
                                      for status in ['PENDING', 'RUNNING']]

class Flow(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    serialization = models.TextField(null=True)
    spec = models.TextField(null=True)
    mission = models.ForeignKey('Mission', null=True, on_delete=models.CASCADE)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in FlowStatuses],
                              default=FlowStatuses.PENDING.name)
    claimed = models.NullBooleanField(null=True, default=False)

    def __str__(self):
        return uuid_model_str(self)

missions_models.append(Flow)

class Job(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    name = models.CharField(null=True, max_length=1024)
    job_spec = JSONField(default=dict)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in JobStatuses],
                              default=JobStatuses.PENDING.name)
    claimed = models.NullBooleanField(null=True, default=False)
    data = models.TextField(null=True)
    error = models.TextField(null=True, blank=True)

missions_models.append(Job)

class FlowJob(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=str_uuid4,
                            editable=False, max_length=64)
    flow = models.ForeignKey(Flow, on_delete=models.CASCADE,
                                 related_name='flow_jobs')
    job = models.ForeignKey('Job', on_delete=models.CASCADE)
    finished = models.NullBooleanField(null=True, default=False)
    meta = JSONField(default=dict)

    def __str__(self):
        return '<{class_name}: {{flow_id: {f_id}, job_id: {j_id}}}>'.format(
            class_name=self.__class__.__name__,
            f_id=self.flow_id,
            j_id=self.job_id)
