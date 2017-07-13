from mc.job_module_utils.job_module_command_dispatcher import (
    JobModuleCommandDispatcher)

def main():
    job_module = generate_job_module()
    def load_job_module(*args, **kwargs): return job_module
    dispatcher = JobModuleCommandDispatcher(load_job_module_fn=load_job_module)
    job = {'key': 'some_job'}
    job_spec = dispatcher.build_jobdir(job=job)
    dispatcher.run_jobdir(jobdir=job_spec['dir'])

def generate_job_module():
    class MyJobModule:
        def build_jobdir(self, *args, **kwargs): print("build_jobdir")
        def run_jobdir(self, *args, **kwargs): print("run_jobdir")
    return MyJobModule()

if __name__ == '__main__': main()
