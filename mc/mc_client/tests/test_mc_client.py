import collections
import inspect
import unittest
from unittest.mock import call, MagicMock

from .. import mc_client


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.dao = MagicMock()
        self.mc_client = mc_client.MissionControlClient(dao=self.dao)

    def generate_mocks(self, n=3):
        return [MagicMock() for i in range(n)]

    def run_basic_dao_dispatch_test(self, fn_name=None, args=None, kwargs=None):
        client_fn = getattr(self.mc_client, fn_name)
        args = args or self.generate_mock_args(fn=client_fn)
        kwargs = kwargs or self.generate_mock_kwargs(fn=client_fn)
        result = client_fn(*args, **kwargs)
        dao_fn = getattr(self.mc_client.dao, fn_name)
        self.assertEqual(dao_fn.call_args, call(*args, **kwargs))
        self.assertEqual(result, dao_fn.return_value)

    def generate_mock_args(self, fn=None):
        argspec_components = self.get_argspec_components(fn=fn)
        args = self.generate_mocks(n=len(argspec_components['arg_names']))
        return args

    def get_argspec_components(self, fn=None):
        argspec = inspect.getargspec(fn)
        num_kwargs = len(argspec.defaults or [])
        arg_names = argspec.args[:(num_kwargs - 1)]
        if num_kwargs: kwarg_names = argspec.args[-num_kwargs:]
        else: kwarg_names = []
        if arg_names and arg_names[0] == 'self':
            arg_names = arg_names[1:]
        return {'arg_names': arg_names, 'kwarg_names': kwarg_names}

    def generate_mock_kwargs(self, fn=None):
        argspec_components = self.get_argspec_components(fn=fn)
        kwargs = {kwarg_name: MagicMock()
                  for kwarg_name in argspec_components['kwarg_names']}
        return kwargs

class CreateQueueTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='create_queue')

class ClaimQueueItemsTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='claim_queue_items')

class CreateFlowTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='create_flow')

class GetFlowsTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='get_flows')

class GetFlowByUUIDTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_client.get_flows = MagicMock(return_value=[MagicMock()])

    def test_wraps_get_flows(self):
        uuid = 'some uuid'
        result = self.mc_client.get_flow_by_uuid(uuid=uuid)
        self.assertEqual(
            self.mc_client.get_flows.call_args,
            call(query={
                'filters': [
                    {'field': 'uuid', 'operator': '=', 'value': uuid}
                ]
            })
        )
        self.assertEqual(result, self.mc_client.get_flows.return_value[0])

class PatchFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_client.patch_flow = MagicMock()
        self.keyed_patches = collections.OrderedDict()
        for i, mock in enumerate(self.generate_mocks()):
            self.keyed_patches[i] = mock

    def _patch_flows(self):
        return self.mc_client.patch_flows(keyed_patches=self.keyed_patches)

    def test_makes_patch_calls(self):
        self._patch_flows()
        self.assertEqual(
            self.mc_client.patch_flow.call_args_list,
            [call(key=key, patches=patches_for_key)
             for key, patches_for_key in self.keyed_patches.items()]
        )

    def test_returns_keyed_results(self):
        result = self._patch_flows()
        self.assertEqual(result,
                         {key: self.mc_client.patch_flow.return_value
                          for key in self.keyed_patches})

class PatchFlowTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='patch_flow')

class CreateJobTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='create_job')

class GetJobsTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='get_jobs')

class GetJobByUUIDTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_client.get_jobs = MagicMock(return_value=[MagicMock()])

    def test_wraps_get_jobs(self):
        uuid = 'some uuid'
        result = self.mc_client.get_job_by_uuid(uuid=uuid)
        self.assertEqual(
            self.mc_client.get_jobs.call_args,
            call(query={
                'filters': [
                    {'field': 'uuid', 'operator': '=', 'value': uuid}
                ]
            })
        )
        self.assertEqual(result, self.mc_client.get_jobs.return_value[0])

class PatchJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_client.patch_job = MagicMock()
        self.keyed_patches = collections.OrderedDict()
        for i, mock in enumerate(self.generate_mocks()):
            self.keyed_patches[i] = mock

    def _patch_jobs(self):
        return self.mc_client.patch_jobs(keyed_patches=self.keyed_patches)

    def test_makes_patch_calls(self):
        self._patch_jobs()
        self.assertEqual(
            self.mc_client.patch_job.call_args_list,
            [call(key=key, patches=patches_for_key)
             for key, patches_for_key in self.keyed_patches.items()]
        )

    def test_returns_keyed_results(self):
        result = self._patch_jobs()
        self.assertEqual(result,
                         {key: self.mc_client.patch_job.return_value
                          for key in self.keyed_patches})

class PatchJobTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='patch_job')

class FlushTestCase(BaseTestCase):
    def test_dispatches_to_dao(self):
        self.run_basic_dao_dispatch_test(fn_name='flush_mc_db')

if __name__ == '__main__':
    unittest.main()


