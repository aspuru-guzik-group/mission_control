import collections
import glob
import json
import os

from mc.a2g2.a2g2_client import a2g2_client as a2g2_client_module


class A2G2_DB_Uploader(object):
    def upload_bulk_files(self, *args, bulk_files_dir=None, cfg=None, **kwargs):
        a2g2_client = self.generate_a2g2_client(cfg=cfg)
        parsed_items = self.parse_bulk_files_dir(
            bulk_files_dir=bulk_files_dir)
        return self.upload_chemthings(
            chemthings=parsed_items.get('chemthings', []),
            a2g2_client=a2g2_client)

    def generate_a2g2_client(self, cfg=None):
        a2g2_client_cfg_json = self.get_cfg_value(
            cfg=cfg, key='a2g2_client_cfg_json', default='{}')
        a2g2_client_cfg = json.loads(a2g2_client_cfg_json)
        a2g2_client = a2g2_client_module.A2G2_Client(**a2g2_client_cfg)
        return a2g2_client

    def get_cfg_value(self, cfg=None, key=None, default=None):
        if key in os.environ: return os.environ[key]
        return (cfg or {}).get(key, default)

    def parse_bulk_files_dir(self, bulk_files_dir=None):
        parsed_items = collections.defaultdict(list)
        bulk_file_paths = self.get_bulk_file_paths(
            bulk_files_dir=bulk_files_dir)
        for bulk_file_path in bulk_file_paths:
            bulk_file_items = self.parse_bulk_file(
                bulk_file_path=bulk_file_path)
            for item_type, items_for_type in bulk_file_items.items():
                parsed_items[item_type].extend(items_for_type)
        return parsed_items

    def get_bulk_file_paths(self, bulk_files_dir=None):
        return glob.glob(bulk_files_dir + '/*.chemthings.bulk')

    def parse_bulk_file(self, bulk_file_path=None):
        bulk_file_items = collections.defaultdict(list)
        with open(bulk_file_path, 'r') as bulk_file:
            for line_num, line in enumerate(bulk_file.readlines()):
                try:
                    item = json.loads(line.strip())
                    # @TODO: for now, we assume everything is a chemthing.
                    # in the future, could make this more general, or allow
                    # for other actions.
                    bulk_file_items['chemthings'].append(item)
                except Exception as error:
                    raise Exception(("Invalid bulk_file line #{line_num}:\n"
                                     "Line: '{line}'\n."
                                     "Error: '{error}'").format(
                                         line_num=line_num,
                                         line=line,
                                         error=error))
        return bulk_file_items

    def upload_chemthings(self, chemthings=None, a2g2_client=None):
        results = []
        for chemthing in chemthings:
            result = a2g2_client.create_chemthing(chemthing=chemthing)
            results.append(result)
        return results

def generate_uploader():
    return A2G2_DB_Uploader()

def load_from_input_dir(*args, input_dir=None, cfg=None, **kwargs):
    uploader = A2G2_DB_Uploader()
    submission = json.load(open(os.path.join(input_dir, 'submission.json')))
    outputs_dir = os.path.join(input_dir, submission['outputs_dir'])
    result = {}
    bulk_files_dirs = os.listdir(outputs_dir)
    for bulk_files_dir in bulk_files_dirs: 
        full_path = os.path.join(outputs_dir, bulk_files_dir)
        result[bulk_files_dir] = uploader.upload_bulk_files(
            bulk_files_dir=full_path, cfg=cfg)
    return result
