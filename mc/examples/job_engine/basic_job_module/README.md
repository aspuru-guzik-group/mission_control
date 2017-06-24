How do you go from a job spec to actual output? There are several steps.
A simple approach would be to write a python function that runs the job. The advantage of this approach would be its simplicity. However, there would be some disadvantages:
1. Blocking: If it took a long time to run the job, it would block the job_runner from running other jobs.
2. Dependencies: If the job required external dependencies, such as 3rd-party software or libraries, the handler function would need to know how to load those dependencies, without polluting the job_runner's environment. Furthermore, we would need to install these dependencies in the same environment that the job_runner runs in.
3. Environment: If the job needed specific computational resources, it would need to be run in an environment with access to those resources. For example, access to networked disk storage, or to a Spark cluster. This means that the job_runner would also need to run in that same environment.

These disadvantages become problematic for many types of projects. For example, a computational chemistry project often needs to run external software that performs expensive calculations on multi-node clusters. Often these calculations take hours to run. So it wouldn't make sense to run them in the context of the job_runner.

So what are some approaches we can take to address these disadvantages?

One approach is to (1) build a submission spec, and then (2) to submit that submission spec to an external service that actually runs the submission. Then (3) the job_runner periodically polls the external service for the status of the job. When (4) the job_runner sees the job is finished, it updates the job in the mission_control database.

This addresses the disadvantages we encountered with our simple approach.
1. Blocking: We now submit and run jobs asynchrously to external resources. So that job_runner doesn't get blocked on individual jobs.
2. Dependencies: The external resources now manage dependencies, so we don't have manage these things in job_runner. For example, a computational cluster provides an optimized version of a computational chemistry package.
2. Environment: The computational cluster is connected to shared disk, so we don't have to worry about connecting to shared disk in the environment our job_runner runs in.

However, we now have a new set of problems:
1. How do we build submission specs?
2. How do we submit these specs?
3. How do we track submissions?

mission_control provides a framework for helping with these issues, called the Job Engine.

In general, what the job engine does is this:
1. It translates a mission_control job into a submission_spec for a given computation target. E.g. for a specific Slurm cluster.
2. It generates an entrypoint script to run on the computation target.
3. It defines standard checkpoint files to indicate job completion and failure.
4. It defines standard log files to expose stderr and stdout.
5. It exposes job artifacts in a standard way.

That is, it acts an adapter between a mission_control job record and a cluster-specific submission.

The JobEngine is organized around a central JobEngine class, that dispatches to job_modules. This means that users can define custom job_modules without having to modify the central job_engine code.

So let's explore this. First let's define basic job_module.
