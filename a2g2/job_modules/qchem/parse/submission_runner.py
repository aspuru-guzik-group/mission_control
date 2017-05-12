from ...a2g2_common import base_transform_submission_runner
from . import qchem_parser
from .. import constants

def run_job_submission(*args, **kwargs):
    return base_transform_submission_runner.run_job_submission(
        *args,
        transform_fn=qchem_parser.parse_qchem_dir,
        output_key=constants.PARSE_OUTPUTS_KEY,
        **kwargs)
