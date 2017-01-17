import io
import json
import shlex
import unittest
from unittest.mock import patch

from django.conf.urls import url, include
from django.test import TestCase, override_settings

from missions.models import Flow

from ..commands.create_flow_record import Command as CreateFlowRecordCommand
from ..commands import base as base_command


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s/' % BASE_PATH, include('missions.urls')),
]

class BaseTestCase(TestCase):
    def setUp(self):
        self.configure_request_client()
        self.request_client_patcher = patch.object(
            base_command, 'request_client', new=self.client)
        self.request_client_patcher.start()

    def configure_request_client(self):
        for method_name in ['patch', 'post']:
            self.patch_client_method_to_use_json(method_name=method_name)

    def patch_client_method_to_use_json(self, method_name=None):
        orig_method = getattr(self.client, method_name)
        def patched_method(path, data=None, **kwargs):
            return orig_method(path, json.dumps(data), 
                               content_type='application/json', **kwargs)
        setattr(self.client, method_name, patched_method)

    def tearDown(self):
        self.request_client_patcher.stop()

@override_settings(ROOT_URLCONF=__name__)
class CreateFlowRecordE2ETestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_creates_flow_record(self):
        flow_spec = {'flow_type': 'some flow type'}
        parsed_output = self.run_create_flow_record_command(flow_spec=flow_spec)
        self.assert_flow_created(flow=parsed_output['flow'],
                                 flow_spec=flow_spec)

    def run_create_flow_record_command(self, flow_spec=None):
        arg_strs = [
            "--flow_spec_json='%s'" % json.dumps(flow_spec),
            "--flow_server_url='/%s/'" % BASE_PATH,
            "--job_server_url=''",
        ]
        args = shlex.split(' '.join(arg_strs))
        stdout = io.StringIO()
        CreateFlowRecordCommand.run(args=args, streams={'stdout': stdout})
        raw_output = stdout.getvalue()
        parsed_output = json.loads(raw_output)
        return parsed_output

    def assert_flow_created(self, flow=None, flow_spec=None):
        flow_record = Flow.objects.get(uuid=flow['uuid'])
        self.assertEqual(flow_record.spec, json.dumps(flow_spec))

@unittest.skip
class RunFlowTestCase(TestCase):
    def test_runs_flow(self):
        self.fail()
