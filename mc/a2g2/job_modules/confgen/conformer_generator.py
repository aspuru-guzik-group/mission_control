from mc.a2g2.conformer_generators.rdkit_conformer_generator\
        .rdkit_conformer_generator import RDKitConformerGenerator

def generate_conformers(**kwargs):
    generator = RDKitConformerGenerator(**kwargs)
    generator.generate_conformers()
