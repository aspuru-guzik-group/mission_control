import enum
import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel

class JobStatuses(enum.Enum):
    Pending = {'label': 'pending'}
    Claimed = {'label': 'claimed'}
    Completed = {'label': 'completed'}

class Job(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4,
                            editable=False, max_length=64)
    name = models.CharField(null=True, max_length=1024)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in JobStatuses],
                              default=JobStatuses.Pending.name)
