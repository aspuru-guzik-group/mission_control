import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import JSONField

from .constants import JobStatuses

class Job(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4,
                            editable=False, max_length=64)
    name = models.CharField(null=True, max_length=1024)
    spec = JSONField(default=dict)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in JobStatuses],
                              default=JobStatuses.PENDING.name)
