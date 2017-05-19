from django.db.models import Q


def claim_queue_items(queue=None, models=None, db_id=None):
    qset = generate_base_queryset_for_queue(queue=queue, models=models,
                                            db_id=db_id)
    items_to_claim = qset.exclude(get_claimed_filter())
    items_to_return = list(items_to_claim)
    claim_queryset_items(queryset=items_to_claim)
    return items_to_return

def serialize_queue_items(queue=None, items=None, serializers=None):
    serializer = getattr(serializers,
                         '%sSerializer' %queue.queue_spec['item_type'])
    return [serializer(item).data for item in items]

def generate_base_queryset_for_queue(queue=None, models=None, db_id=None):
    ItemModel = getattr(models, queue.queue_spec['item_type'])
    if db_id: queryset = ItemModel.objects.using(db_id).filter()
    else: queryset = ItemModel.objects.filter()
    return queryset

def exclude_claimed_items_from_queryset(queryset=None):
    if hasattr(queryset.model, 'claimed'):
        queryset = queryset.exclude(claimed=True)
    return queryset

def claim_queryset_items(queryset=None):
    if hasattr(queryset.model, 'claimed'): queryset.update(claimed=True)

def get_claimed_filter():
    return Q(claimed=True)
