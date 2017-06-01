import json

from mc.a2g2.a2g2_client import a2g2_client as a2g2_client_module
from ...job_modules import utils as job_module_utils

def generate_a2g2_client(cfg=None):
    a2g2_client_cfg_json = job_module_utils.get_cfg_value(
        cfg=cfg, key='a2g2_client_cfg_json', default='{}')
    a2g2_client_cfg = json.loads(a2g2_client_cfg_json)
    a2g2_client = a2g2_client_module.A2G2_Client(**a2g2_client_cfg)
    return a2g2_client
