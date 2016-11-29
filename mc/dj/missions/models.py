import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel


class Mission(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)

class Task(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)
    mission = models.ForeignKey('Mission', null=True, on_delete=models.CASCADE)

class Workflow(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    runner_key = models.CharField(null=True, max_length=1024)
    mission = models.ForeignKey('Mission', null=True, on_delete=models.CASCADE)
    finished = models.NullBooleanField(null=True)

    @property
    def last_finished_job(self):
        if not self.jobs.exists(): return None
        return self.jobs.filter(finished=True).order_by('-modified')[0]

class WorkflowRunner(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)
