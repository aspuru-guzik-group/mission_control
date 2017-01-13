import io
import json
import re
from django.conf.urls import url, include
from django.test import TestCase, override_settings

from missions.models import Flow

from ..commands.create_flow_record import Command as CreateFlowRecordCommand


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('missions.urls')),
]

class BaseTestCase(TestCase):
    def setUp(self):
        self.configure_request_client()

    def configure_request_client(self):
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch

@override_settings(ROOT_URLCONF=__name__)
class CreateFlowRecordE2ETestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_creates_flow_record(self):
        flow_spec = {'flow_type': 'some flow type'}
        parsed_output = self.run_create_flow_record_command(flow_spec=flow_spec)
        self.assert_flow_created(flow_uuid=parsed_output['flow_uuid'],
                                 flow_spec=flow_spec)

    def run_create_flow_record_command(self, flow_spec=None):
        args = [
            "--flow_spec_json='%s'" % json.dumps(flow_spec),
            "--flow_server_url='%s'" % ("/%s/flows/" % BASE_PATH),
            "--job_server_url=''",
        ]
        stdout = io.StringIO()
        CreateFlowRecordCommand.run(args=args, streams={})
        parsed_output = self.parse_command_output(output=stdout.getvalue())
        return parsed_output

    def parse_command_output(self, output=None):
        match = re.search(r"Created flow with uuid '(.+)'", output)
        if not match: parsed_output = None
        else: parsed_output = {'uuid': match.group(1)}
        return parsed_output

    def assert_flow_created(self, flow_uuid=None, flow_spec=None):
        flow_record = Flow.objects.get(uuid=flow_uuid)
        self.assertEqual(flow_record.spec, json.dumps(flow_spec))

class RunFlowTestCase(TestCase):
    def test_runs_flow(self):
        self.fail()
