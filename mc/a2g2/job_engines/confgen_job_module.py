from mc.a2g2.a2g2_utils.conformer_generators.rdkit_conformer_generator\
        .rdkit_conformer_generator import RDKitConformerGenerator

def execute_job(*args, job=None, cfg=None, output_dir=None, **kwargs):
    command = job['job_spec']['command']
    if command == 'confgen':
        confgen_kwargs = {**job['job_spec']['kwargs'], 'output_dir': output_dir}
        generator = RDKitConformerGenerator(**confgen_kwargs)
        generator.generate_conformers()
    elif command == 'parse':
        pass
    elif command == 'load':
        pass
