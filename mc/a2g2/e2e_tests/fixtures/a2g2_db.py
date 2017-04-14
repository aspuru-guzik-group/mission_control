import json
import sys


class A2G2_DB_Fixtures(object):
    def __init__(self):
        pass

    def run_fake_query(self):
        fake_result = {
            'query': {},
            'hits': [1, 2, 3]
        }
        sys.stdout.write(json.dumps(fake_result))
