import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel


class Job(TimeStampedModel):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4,
                            editable=False)
    name = models.CharField(null=True, max_length=1024)
