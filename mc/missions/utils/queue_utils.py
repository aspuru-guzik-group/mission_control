from django.apps import apps as _apps

from ..serializers import missions_serializers


def get_serialized_queue_items(queue=None, query_params=None):
    return serialize_queue_items(
        queue=queue,
        queue_items=get_queue_items(queue=queue, query_params=query_params))

def serialize_queue_items(queue=None, queue_items=None):
    serializer = get_serializer_for_queue(queue=queue)
    return [serializer(queue_item).data for queue_item in queue_items]

def get_serializer_for_queue(queue=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue.queue_spec)
    return missions_serializers[RootModel]

def get_root_model_for_queue_spec(queue_spec=None):
    return _apps.get_model(queue_spec['root_model_spec'])

def get_queue_items(queue=None, query_params=None):
    base_queryset = generate_base_queryset_for_queue_spec(
        queue_spec=queue.queue_spec)
    queue_items = base_queryset
    return queue_items

def generate_base_queryset_for_queue_spec(queue_spec=None):
    RootModel = get_root_model_for_queue_spec(queue_spec=queue_spec)
    return RootModel.objects.filter()


