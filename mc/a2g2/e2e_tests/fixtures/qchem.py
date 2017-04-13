import argparse
import gzip
import os

from mc.a2g2.job_modules.qchem.parse import tests as qchem_parse_tests


class QChemFixtures(object):
    def __init__(self):
        self.molecule = {
            'charge': 0,
            'multiplicity': 1,
            'atoms': [
                {'element': 'C', 'x': i, 'y': i, 'z': i}
                for i in range(3)
            ]
        }
        self.precursors = {'precursor1': True}

    def generate_qchem_output_file_content(self):
        example_confgen_output_path = os.path.join(
            os.path.dirname(qchem_parse_tests.__file__), 'example_qchem.out.gz')
        return gzip.open(example_confgen_output_path, 'rb').read()
    
    @classmethod
    def handle_qchem_command(cls, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument('input_file')
        parser.add_argument('output_file')
        parsed_args, extra_args = parser.parse_known_args(args=args)
        fixtures = QChemFixtures()
        open(parsed_args.output_file,'wb').write(
            fixtures.generate_qchem_output_file_content())
