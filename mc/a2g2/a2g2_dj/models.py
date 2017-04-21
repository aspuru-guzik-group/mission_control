import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.contrib.postgres.fields import JSONField

from taggit.managers import TaggableManager


a2g2_dj_models = []

class ChemThing(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4, max_length=64)
    keys = JSONField(default=dict, null=True)
    types = JSONField(default=dict, null=True)
    props = JSONField(default=dict, null=True)
    precursors = JSONField(default=dict, null=True)
    ancestors = JSONField(default=dict, null=True)
    tags = TaggableManager()

a2g2_dj_models.append(ChemThing)
