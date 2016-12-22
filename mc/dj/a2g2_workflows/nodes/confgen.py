from .job_wrapper import JobWrapperNode


class ConfgenNode(JobWrapperNode):
    def __init__(self, *args, data=None, **kwargs):
        self.validate_data(data=data)
        super().__init__(*args, data=data, **kwargs)

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
