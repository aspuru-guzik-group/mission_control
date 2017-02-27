import os
import textwrap
import sys

from .base_command import BaseCommand
from . import chemthing_utils


def execute_job(*args, job=None, cfg=None, output_dir=None, ctx_dir=None,
                **kwargs):
    Command().handle(*args, job=job, cfg=cfg, output_dir=None, ctx_dir=ctx_dir,
                     **kwargs)

class QChemJobEngine(object):
    def execute_job(self, *args, job=None, cfg=None, output_dir=None,
                    ctx_dir=None, **kwargs):
        qchem_dir_meta = self.generate_qchem_dir(job=job, output_dir=output_dir)
        self.execute_dir(dir_meta=qchem_dir_meta)

    def generate_qchem_dir(self, job=None, output_dir=None):
        chemthing = job['data']['input']['chemthing']

        def _write_dir_file(filename=None, content=None):
            with open(os.path.join(output_dir, filename), 'w') as f:
                f.write(content)

        input_file_name = 'qchem.inp'

        _write_dir_file(
            filename=input_file_name,
            content=self._generate_qchem_input_content(job=job)
        )
        _write_dir_file(
            filename='coords.xyz',
            content=chemthing_utils.from_cml(cml=chemthing['cml'],
                                             output_format='xyz')
        )
        entrypoint_file_name = 'entrypoint.sh'
        _write_dir_file(
            filename=entrypoint_file_name,
            content=textwrap.dedent(
                '''
                #!/bin/bash
                INPUT_FILE="{input_file_name}"
                OUTPUT_FILE="{output_file_name}"
                qchem $INPUT_FILE $OUTPUT_FILE
                '''
            ).format(
                input_file_name=input_file_name,
                output_file_name='qchem.out'
            )
        )
        dir_meta = {
            'dir': output_dir,
            'entrypoint': entrypoint_file_name
        }
        return dir_meta

    def _generate_qchem_input_content(self, job=None):
        return textwrap.dedent(
            '''
            $molecule
            READ coords.xyz
            $end

            $rem
            MEM_TOTAL 12000
            MEM_STATIC 2500
            JOBTYPE OPT
            GEOM_OPT_MAX_CYCLES 200
            GEOM_OPT_PRINT 2

            SYMMETRY FALSE

            !BP86
            EXCHANGE B
            CORRELATION P86

            BASIS 6-31G*
            PURECART 1111

            SCF_GUESS CORE
            XC_GRID 1
            MAX_SCF_CYCLES 200

            PRINT_ORBITALS FALSE
            SCF_PRINT 1
            SCF_FINAL_PRINT 1
            POP_MULLIKEN 0

            $end
            '''
        )

    def test_returns_expected_dir_meta(self):
        dir_meta = self.engine.generate_qchem_dir(job=self.job,
                                                  output_dir=self.output_dir)
        expected_dir_meta = {
            'dir': self.output_dir,
            'entrypoint': 'entrypoint.sh',
        }
        self.assertEqual(dir_meta, expected_dir_meta)

    def execute_dir(self, dir_meta=None):
        pass

class Command(BaseCommand):
    help = 'qchem_job'

    def execute_job(self, *args, job=None, cfg=None, ctx_dir=None, 
                    output_dir=None, **kwargs):
        engine = QChemJobEngine()
        engine.execute_job(job=job, cfg=cfg, ctx_dir=ctx_dir,
                           output_dir=output_dir)

if __name__ == '__main__':
    command = Command().execute(args=sys.argv[1:])
