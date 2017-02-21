import argparse
import json
import sys

from .rdkit_conformer_generator import RDKitConformerGenerator

class Command(object):
    help = 'Generates conformers using RDKit and OpenBabel'

    def setup_streams(self, stdin=sys.stdin, stdout=sys.stdout,
                      stderr=sys.stderr):
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

    def add_arguments(self, parser):
        parser.add_argument(
            '--params_file',
            help=("Input parameter file, in .json format. Parameters from"
                  " this file will override parameters provided on the command" 
                  " line."),
            type=str
        )
        parser.add_argument(
            '--output_dir',
            type=str
        )
        parser.add_argument(
            '--forcefield_id',
            help="id of forcefield model to use.",
            choices=RDKitConformerGenerator.forcefields,
            default='mmff'
        )
        parser.add_argument(
            '--smiles',
            type=str
        )
        parser.add_argument(
            '--output_limit',
            help=("Maximum number of conformers to output. Cut-off is"
                  " performed after any other filtering."),
            type=int,
            default=20
        )
        parser.add_argument(
            '--candidate_limit',
            help=("Number of initial candidate conformers to generate. These"
                  " candidates will be filtered and clustered to produce the"
                  " final output set. Setting this number higher, in"
                  " combination with the output_limit can potentially product a"
                  " greater number of conformers, but at the expense of longer"
                  " compute times."),
            type=int,
            default=200
        )
        parser.add_argument(
            '--prune_rms_thresh',
            help=("RDKit's 'pruneRmsThresh' parameter. See: "
                  "http://www.rdkit.org/Python_Docs/rdkit.Chem.rdDistGeom-module.html."
                  ),
            type=float,
            default=0.01
        )
        parser.add_argument(
            '--cluster_energy_range',
            help=("The maximum difference between"
                  " (A) a cluster's minimum energy, and"
                  " (B) the minimum of energy of all canddiate conformers."),
            type=float,
            default=5.0
        )
        parser.add_argument(
            '--cluster_energy_radius',
            help=("The maximum difference between"
                  " (A) <a cluster's minimum energy>, and"
                  " (B) <the energy of candidate conformer under consideration"
                  " for the cluster>."),
            type=float,
            default=5.0
        )
        parser.add_argument(
            '--cluster_energy_auto_accept_radius',
            help=("If the difference between"
                  " (A) <a cluster's minimum energy>, and"
                  " (B) <the energy of candidate conformer under consideration"
                  " for the cluster> is less than this value, it will be"
                  " automatically added to the cluster. This is intended to"
                  " prioritize energy similarity when clustering, even if rms"
                  " values differ."),
            type=float,
            default=1e-3
        )
        parser.add_argument(
            '--cluster_rms_radius',
            help=("The maximum value of RDKit's 'GetBestRMS' function,"
                  " computed for"
                  " (A) <a cluster's seed>, and"
                  " (B) <a candidate conformer under consideration"
                  " for the cluster>."),
            type=float,
            default=0.1
        )
        parser.add_argument(
            '--fallback_to_align_for_rms',
            help=("Whether to use OpenBabel's 'align' method if rms cannot be"
                  " determined via OpenBabel's 'obfit' method."),
            action='store_true'
        )
        parser.add_argument(
            '--xyz_filename_tpl',
            help=("A python string format template to use for creating"
                  " conformer file names."),
            type=str,
            default='conf_{index}.xyz',
        )

    def handle(self, *args, **options):
        params = options
        params_file = params.get('params_file', None)
        if params_file: params.update(json.load(open(params_file)))
        self.generator = RDKitConformerGenerator(**params)
        self.generator.generate_conformers()

if __name__ == '__main__':
    command = Command()
    command.setup_streams()
    parser = argparse.ArgumentParser(description=command.help)
    command.add_arguments(parser)
    args = vars(parser.parse_args())
    command.handle(**args)
