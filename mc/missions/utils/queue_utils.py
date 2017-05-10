from django.db.models import Q


def claim_queue_items(queue=None, models=None):
    qset = generate_base_queryset_for_queue(queue=queue, models=models)
    items_to_claim = qset.exclude(get_claimed_filter())
    items_to_return = list(items_to_claim)
    claim_queryset_items(queryset=items_to_claim)
    return items_to_return

def serialize_queue_items(queue=None, items=None, serializers=None):
    serializer = serializers[queue.queue_spec['item_type']]
    return [serializer(item).data for item in items]

def generate_base_queryset_for_queue(queue=None, models=None):
    ItemModel = models[queue.queue_spec['item_type']]
    return ItemModel.objects.filter()

def exclude_claimed_items_from_queryset(queryset=None):
    if hasattr(queryset.model, 'claimed'):
        queryset = queryset.exclude(claimed=True)
    return queryset

def claim_queryset_items(queryset=None):
    if hasattr(queryset.model, 'claimed'): queryset.update(claimed=True)

def get_claimed_filter():
    return Q(claimed=True)
