import io
import json
import shlex

from django.test import override_settings

from missions.models import Flow
from ..commands.create_flow_record import Command as CreateFlowRecordCommand
from . import e2e_utils

urlpatterns = e2e_utils.urlpatterns
assert urlpatterns

@override_settings(ROOT_URLCONF=__name__)
class CreateFlowRecordE2ETestCase(e2e_utils.BaseTestCase):
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
            "--mc_server_url='/%s/'" % e2e_utils.BASE_PATH,
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
