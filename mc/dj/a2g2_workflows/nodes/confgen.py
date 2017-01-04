from .job_wrapper import JobWrapperNode


class Confgen_Node(JobWrapperNode):
    def __init__(self, *args, **kwargs): super().__init__(self, *args, **kwargs)

    def validate_data(self, data=None):
        assert data['input']['smiles'] is not None
        assert data['input']['confgen_spec'] is not None

    def get_job_input(self):
        return {
            'job_type': 'confgen',
            'job_spec': {
                'smiles': self.data['input']['smiles'],
                'confgen': self.data['input']['confgen_spec']
            }
        }

class Confgen_Parse_Node(JobWrapperNode):
    def __init__(self, *args, **kwargs): super().__init__(self, *args, **kwargs)

    def validate_data(self, data=None):
        assert data['input']['parse']['dir'] is not None

    def get_job_input(self):
        return {
            'job_type': 'confgen_parse',
            'job_spec': {
                'confgen_dir': self.data['input']['parse']['dir']
            }
        }
