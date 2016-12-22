from .job_wrapper import JobWrapperNode


class B3LYPNode(JobWrapperNode):
    def __init__(self, *args, data=None, **kwargs):
        self.validate_data(data=data)
        super().__init__(*args, data=data, **kwargs)

    def validate_data(self, data=None):
        assert data['input']['xyz'] is not None
        assert data['input']['b3lyp_spec'] is not None

    def get_job_input(self):
        return {
            'job_type': 'b3lyp',
            'job_spec': {
                'xyz': self.data['input']['xyz'],
                'b3lyp': self.data['input']['b3lyp_spec']
            }
        }
