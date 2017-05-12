import logging

from mc.missions.utils import queue_utils as _queue_utils


class DjDao(object):
    def __init__(self, models=None, serializers=None, queue_utils=None,
                 logger=None):
        self.models = models
        self.serializers = serializers
        self.queue_utils = queue_utils or _queue_utils
        self.logger = logger or logging

    def create_flow(self, flow_kwargs=None):
        flow_model = self.models['Flow'].objects.create(**flow_kwargs)
        return self.serialize_flow_model(flow_model=flow_model)

    def serialize_flow_model(self, flow_model=None):
        return self.serializers['Flow'](flow_model).data

    def get_flows(self, query=None):
        flow_models = self.get_flow_models(query=query)
        return [self.serialize_flow_model(flow_model=flow_model)
                for flow_model in flow_models]

    def get_flow_models(self, query=None):
        flow_models = []
        queryset = self.flow_query_to_queryset(query=query)
        if queryset: flow_models = queryset.all()
        return flow_models
 
    def flow_query_to_queryset(self, query=None):
        if query: self.validate_query(query=query)
        qs = self.models['Flow'].objects
        return qs.filter(**self.get_dj_filter_kwargs_for_query(query=query))

    def validate_query(self, query=None):
        for _filter in query.get('filters' or []):
            self.validate_query_filter(_filter=_filter)

    def validate_query_filter(self, _filter=None):
        valid_operators = ['=', 'IN']
        if _filter['operator'] not in valid_operators:
            raise Exception(("Invalid operator '{operator}',"
                             " valid operators are: {valid_operators}").format(
                                 operator=_filter['operator'],
                                 valid_operators=valid_operators))

    def get_dj_filter_kwargs_for_query(self, query=None):
        filter_kwargs = {}
        if not query: return filter_kwargs
        for _filter in query.get('filters', []):
            filter_kwargs.update(
                self.query_filter_to_dj_filter_kwargs(_filter=_filter))
        return filter_kwargs

    def query_filter_to_dj_filter_kwargs(self, _filter=None):
        dj_filter_kwargs = {}
        kwarg_name = _filter['field']
        if _filter['operator'] != '=':
            kwarg_name += _filter['operator'].lower()
        dj_filter_kwargs[kwarg_name] = _filter['value']
        return dj_filter_kwargs

    def patch_flow(self, key=None, patches=None):
        flow_model = self.models['Flow'].objects.get(uuid=key)
        self.patch_model(model=flow_model, patches=patches)
        return self.serialize_flow_model(flow_model=flow_model)

    def patch_model(self, model=None, patches=None):
        patches = patches or {}
        for k, v in patches.items():
            setattr(model, k, v)
        model.save()

    def create_job(self, job_kwargs=None):
        job_model = self.models['Job'].objects.create(**job_kwargs)
        return self.serialize_job_model(job_model=job_model)

    def serialize_job_model(self, job_model=None):
        return self.serializers['Job'](job_model).data

    def get_jobs(self, query=None):
        job_models = self.get_job_models(query=query)
        return [self.serialize_job_model(job_model=job_model)
                for job_model in job_models]

    def get_job_models(self, query=None):
        job_models = []
        queryset = self.job_query_to_queryset(query=query)
        if queryset: job_models = queryset.all()
        return job_models
 
    def job_query_to_queryset(self, query=None):
        if query: self.validate_query(query=query)
        qs = self.models['Job'].objects
        return qs.filter(**self.get_dj_filter_kwargs_for_query(query=query))

    def patch_job(self, key=None, patches=None):
        job_model = self.models['Job'].objects.get(uuid=key)
        self.patch_model(model=job_model, patches=patches)
        return self.serialize_job_model(job_model=job_model)

    def create_queue(self, queue_kwargs=None):
        queue_model = self.models['Queue'].objects.create(**queue_kwargs)
        return self.serialize_queue_model(queue_model=queue_model)

    def serialize_queue_model(self, queue_model=None):
        return self.serializers['Queue'](queue_model).data

    def flush_mc_db(self):
        for model_name in ['Job', 'Flow', 'Queue']:
            self.models[model_name].objects.all().delete()

    def claim_queue_items(self, queue_key=None, params=None):
        queue_model = self.models['Queue'].objects.get(uuid=queue_key)
        claimed_item_models = self.queue_utils.claim_queue_items(
            queue=queue_model, models=self.models)
        return {
            'items': self.queue_utils.serialize_queue_items(
                queue=queue_model, items=claimed_item_models,
                serializers=self.serializers)
        }
