import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import ArrayField, JSONField


a2g2_dj_models = []

class ChemThing(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4,
                            editable=False, max_length=64)
    types = ArrayField(models.CharField(max_length=128, blank=True),
                       default=list)
    props = JSONField(null=True)

a2g2_dj_models.append(ChemThing)
