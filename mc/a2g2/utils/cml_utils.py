import pybel


def from_cml(cml=None, output_format=None):
    return pybel.readstring('cml', cml).write(output_format)

def to_cml(input=None, input_format=None):
    return pybel.readstring(input_format, input).write('cml')
