import contextlib
import importlib
import logging


class DjDao(object):
    DEFAULT_DB_CFG = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': ':memory:'
        }
    }
    MC_APP = 'mc.missions'
    INSTALLED_APPS_CFG = [MC_APP]

    def __init__(self, db_cfg=None, mc_modules=None, logger=None):
        self.logger = logger or logging
        self.db_cfg = db_cfg or self.DEFAULT_DB_CFG
        self._mc_modules = mc_modules

    @contextlib.contextmanager
    def get_mc_modules_ctx(self):
        """This is a context manager that sets self.mc_modules to either
        the mc_modules provided to __init__, or modules loaded from
        the mc django app, in a django context"""
        if hasattr(self, 'mc_modules'): yield
        elif self._mc_modules is not None:
            self.mc_modules = self._mc_modules
            yield
        else:
            with self.get_dj_ctx():
                self.mc_modules = {
                    module_name: self.import_mc_module(module_name)
                    for module_name in ['models', 'serializers']
                }
                self.mc_modules['queue_utils'] = self.import_mc_module(
                    'utils.queue_utils')
                yield
        if hasattr(self, 'mc_modules'): del self.mc_modules

    def import_mc_module(self, module_name=None):
        prefixed_module_name = '{mc_app}.{module_name}'.format(
            mc_app=self.MC_APP, module_name=module_name)
        return importlib.import_module(prefixed_module_name)

    def get_dj_ctx(self):
        """ This is some hackiness to setup a temporary django environment.
        It patches django internals just enough to call django.setup().
        """
        from unittest.mock import patch
        import sys
        from django.apps import registry as _dj_apps_registry
        from django.conf import LazySettings
        with patch.dict(sys.modules, {_dj_apps_registry.__name__: None}):
            _apps = _dj_apps_registry.Apps(installed_apps=None)
        _settings = LazySettings()
        exit_stack = contextlib.ExitStack()
        for ctx_mgr in [
            patch('django.conf.settings', new=_settings),
            patch('django.apps.apps', new=_apps)
        ]: exit_stack.enter_context(ctx_mgr)
        from django.conf import settings
        settings.configure(
            DATABASES=self.db_cfg,
            INSTALLED_APPS=self.INSTALLED_APPS_CFG
        )
        import django
        django.setup(set_prefix=False)
        return exit_stack

    def create_flow(self, flow_kwargs=None):
        with self.get_mc_modules_ctx():
            flow_model = self.mc_modules['models'].Flow.objects.create(
                **flow_kwargs)
            return self.serialize_flow_model(flow_model=flow_model)

    def serialize_flow_model(self, flow_model=None):
        with self.get_mc_modules_ctx():
            return self.mc_modules['serializers'].FlowSerializer(
                flow_model).data

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
        with self.get_mc_modules_ctx():
            if query: self.validate_query(query=query)
            qs = self.mc_modules['models'].Flow.objects
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
        with self.get_mc_modules_ctx():
            flow_model = self.mc_modules['models'].Flow.objects.get(uuid=key)
            self.patch_model(model=flow_model, patches=patches)
            return self.serialize_flow_model(flow_model=flow_model)

    def patch_model(self, model=None, patches=None):
        patches = patches or {}
        for k, v in patches.items(): setattr(model, k, v)
        model.save()

    def create_job(self, job_kwargs=None):
        with self.get_mc_modules_ctx():
            job_model = self.mc_modules['models'].Job.objects.create(
                **job_kwargs)
            return self.serialize_job_model(job_model=job_model)

    def serialize_job_model(self, job_model=None):
        with self.get_mc_modules_ctx():
            return self.mc_modules['serializers'].JobSerializer(job_model).data

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
        with self.get_mc_modules_ctx():
            if query: self.validate_query(query=query)
            qs = self.mc_modules['models'].Job.objects
            return qs.filter(**self.get_dj_filter_kwargs_for_query(query=query))

    def patch_job(self, key=None, patches=None):
        with self.get_mc_modules_ctx():
            job_model = self.mc_modules['models'].Job.objects.get(uuid=key)
            self.patch_model(model=job_model, patches=patches)
            return self.serialize_job_model(job_model=job_model)

    def create_queue(self, queue_kwargs=None):
        with self.get_mc_modules_ctx():
            queue_model = self.mc_modules['models'].Queue.objects.create(
                **queue_kwargs)
            return self.serialize_queue_model(queue_model=queue_model)

    def serialize_queue_model(self, queue_model=None):
        with self.get_mc_modules_ctx():
            return self.mc_modules['serializers'].QueueSerializer(
                queue_model).data

    def claim_queue_items(self, queue_key=None, params=None):
        with self.get_mc_modules_ctx():
            queue_model = self.mc_modules['models'].Queue.objects.get(
                uuid=queue_key)
            queue_utils = self.mc_modules['queue_utils']
            claimed_item_models = queue_utils.claim_queue_items(
                queue=queue_model, models=self.mc_modules['models'])
            return {
                'items': queue_utils.serialize_queue_items(
                    queue=queue_model, items=claimed_item_models,
                    serializers=self.mc_modules['serializers'])
            }

    def flush_mc_db(self):
        with self.get_mc_modules_ctx():
            for model_name in ['Job', 'Flow', 'Queue']:
                ModelCls = getattr(self.mc_modules['models'], model_name)
                ModelCls.objects.all().delete()
