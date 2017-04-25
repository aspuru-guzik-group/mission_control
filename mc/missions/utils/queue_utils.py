from django.apps import apps as _apps

from ..serializers import missions_serializers


def claim_queue_items(queue=None, query_params=None):
    qset = generate_base_queryset_for_queue_spec(queue_spec=queue.queue_spec)
    update_kwargs = generate_claim_update_kwargs_for_queue(queue=queue)
    if update_kwargs: qset.update(**update_kwargs)
    return qset

def serialize_queue_items(queue=None, queue_items=None):
    serializer = get_serializer_for_queue(queue=queue)
    return [serializer(queue_item).data for queue_item in queue_items]

def get_serializer_for_queue(queue=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue.queue_spec)
    return missions_serializers[RootModel]

def get_root_model_for_queue_spec(queue_spec=None):
    return _apps.get_model(queue_spec['root_model_spec'])

def generate_claim_update_kwargs_for_queue(queue=None):
    update_kwargs = {}
    RootModel = get_root_model_for_queue_spec(queue_spec=queue.queue_spec)
    if hasattr(RootModel, 'claimed'): update_kwargs['claimed'] = True
    return update_kwargs

def generate_base_queryset_for_queue_spec(queue_spec=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue_spec)
    return RootModel.objects.filter()


