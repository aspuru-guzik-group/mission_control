from a2g2.a2g2_utils.conformer_generators.rdkit_conformer_generator\
        .rdkit_conformer_generator import RDKitConformerGenerator

def execute_job(*args, job=None, cfg=None, output_dir=None, **kwargs):
    params = {**job['job_spec']['confgen'], 'output_dir': output_dir}
    generator = RDKitConformerGenerator(**params)
    generator.generate_conformers()
