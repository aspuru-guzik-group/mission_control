import enum
import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import JSONField

def uuid_model_str(instance):
    return '<{class_name}: {uuid}>'.format(
        class_name=instance.__class__.__name__,
        uuid=instance.uuid)

class Mission(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)

    def __str__(self):
        return uuid_model_str(self)

class WorkflowStatuses(enum.Enum):
    Pending = {'label': 'pending'}
    Running = {'label': 'running'}
    Completed = {'label': 'completed'}

class Workflow(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    spec = models.ForeignKey('WorkflowSpec', null=True,
                             on_delete=models.CASCADE)
    serialization = models.TextField(null=True)
    mission = models.ForeignKey('Mission', null=True, on_delete=models.CASCADE)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in WorkflowStatuses],
                              default=WorkflowStatuses.Pending.name)

    def __str__(self):
        return uuid_model_str(self)

class WorkflowJob(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    workflow = models.ForeignKey(Workflow, on_delete=models.CASCADE,
                                 related_name='workflow_jobs')
    job = models.ForeignKey('jobs.Job', on_delete=models.CASCADE)
    finished = models.NullBooleanField(null=True, default=False)
    meta = JSONField(default=dict)

    def __str__(self):
        return '<{class_name}: {{workflow_id: {wf}, job_id: {j}}}>'.format(
            class_name=self.__class__.__name__,
            wf=self.workflow_id,
            j=self.job_id)

class WorkflowSpec(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    key = models.CharField(null=True, max_length=1024)
    path = models.CharField(null=True, max_length=1024)
    label = models.CharField(null=True, max_length=1024)
    error = models.CharField(null=True, max_length=1024)
    serialization = models.TextField(null=True)

    def __str__(self):
        return (
            '<{class_name}: {{label: {label}, key: {key}, uuid: {uuid}}}>'
        ).format(
            class_name=self.__class__.__name__,
            label=self.label,
            key=self.key,
            uuid=self.uuid)