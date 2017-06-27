from mc.job_engines.job_engine import JobEngine

def main():
    job_module = setup_job_module()
    job_module_loader = setup_job_module_loader(job_module=job_module)
    job_engine = JobEngine(job_module_loader=job_module_loader)
    job = {'key': 'some_job'}
    submission_meta = job_engine.build_job_submission(job=job)
    job_engine.run_job_submission(submission_dir=submission_meta['dir'])

def setup_job_module():
    class MyJobModule:
        def build_job_submission(self, *args, **kwargs): 
            print("build_job_submission")
        def run_job_submission(self, *args, **kwargs):
            print("run_job_submission")
    return MyJobModule()

def setup_job_module_loader(job_module=None):
    class MyJobModuleLoader:
        def load_job_module(*args, **kwargs): return job_module
    return MyJobModuleLoader()

if __name__ == '__main__': main()
