from mc.job_module_utils.dispatcher import JobModuleCommandDispatcher

def main():
    job_module = setup_job_module()
    job_module_loader = setup_job_module_loader(job_module=job_module)
    dispatcher = JobModuleCommandDispatcher(job_module_loader=job_module_loader)
    job = {'key': 'some_job'}
    jobdir_meta = dispatcher.build_jobdir(job=job)
    dispatcher.run_jobdir(jobdir=jobdir_meta['dir'])

def setup_job_module():
    class MyJobModule:
        def build_jobdir(self, *args, **kwargs): print("build_jobdir")
        def run_jobdir(self, *args, **kwargs): print("run_jobdir")
    return MyJobModule()

def setup_job_module_loader(job_module=None):
    class MyJobModuleLoader:
        def load_job_module(*args, **kwargs): return job_module
    return MyJobModuleLoader()

if __name__ == '__main__': main()
