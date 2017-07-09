Now we try loading from an external job module. This is more likely to be the scenario that users encounter.

By convention, job engine's default job_module_loader gets the job module name from job['job_type'].
It then tries to load the module from that path.
