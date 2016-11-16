import enum
import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel


class Job(TimeStampedModel):
    class STATUSES(enum.Enum):
        PENDING = {'label': 'pending'}
        CLAIMED = {'label': 'claimed'}

    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)
    status = models.CharField(null=True, max_length=32,
                              choices=[(status.name, status.value['label'])
                                       for status in STATUSES],
                              default=STATUSES.PENDING.name)
