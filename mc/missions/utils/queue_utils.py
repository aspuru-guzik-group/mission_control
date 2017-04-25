from django.apps import apps as _apps
from django.db.models import Q

from ..serializers import missions_serializers


def claim_queue_items(queue=None, query_params=None):
    qset = generate_base_queryset_for_queue_spec(queue_spec=queue.queue_spec)
    to_claim = qset.exclude(get_claimed_filter())
    to_return = list(to_claim)
    claim_queryset_items(queryset=to_claim)
    return to_return

def serialize_queue_items(queue=None, queue_items=None):
    serializer = get_serializer_for_queue(queue=queue)
    return [serializer(queue_item).data for queue_item in queue_items]

def get_serializer_for_queue(queue=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue.queue_spec)
    return missions_serializers[RootModel]

def get_root_model_for_queue_spec(queue_spec=None):
    return _apps.get_model(queue_spec['root_model_spec'])

def generate_base_queryset_for_queue_spec(queue_spec=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue_spec)
    return RootModel.objects.filter()

def exclude_claimed_items_from_queryset(queryset=None):
    if hasattr(queryset.model, 'claimed'):
        queryset = queryset.exclude(claimed=True)
    return queryset

def claim_queryset_items(queryset=None):
    if hasattr(queryset.model, 'claimed'): queryset.update(claimed=True)

def get_claimed_filter():
    return Q(claimed=True)
