import uuid
from django.db import models
from django_extensions.db.models import TimeStampedModel

from jsonfield import JSONField

import taggit.managers
import taggit.models


a2g2_dj_models = []

class TaggedChemThing(taggit.models.TaggedItemBase):
    content_object = models.ForeignKey('ChemThing')

class ChemThing(TimeStampedModel):
    uuid = models.CharField(primary_key=True, default=uuid.uuid4, max_length=64)
    keys = JSONField(default=dict, null=True)
    types = JSONField(default=dict, null=True)
    props = JSONField(default=dict, null=True)
    precursors = JSONField(default=dict, null=True)
    ancestors = JSONField(default=dict, null=True)
    tags = taggit.managers.TaggableManager(through=TaggedChemThing)

a2g2_dj_models.append(ChemThing)
