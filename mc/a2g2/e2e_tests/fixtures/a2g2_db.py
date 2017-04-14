import json
import sys


class A2G2_DB_Fixtures(object):
    def __init__(self):
        pass

    def run_fake_conformer_query(self, num_results=3):
        fake_result = {
            'query': {},
            'hits': [
                {
                    'types': {
                        'a2g2:type:mol3d': True
                    },
                    'props': {
                        'a2g2:prop:atoms': self.generate_fake_atoms()
                    }
                } for i in range(num_results)
            ]
        }
        sys.stdout.write(json.dumps(fake_result))

    def generate_fake_atoms(self):
        return [{'element': 'c', 'x': i, 'y': i, 'z': i} for i in range(3)]
