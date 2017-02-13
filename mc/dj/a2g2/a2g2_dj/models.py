import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import JSONField


a2g2_dj_models = []

class ChemThing(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4,
                            editable=False, max_length=64)
    cml = models.TextField(null=True)
    props = JSONField(null=True)

a2g2_dj_models.append(ChemThing)
