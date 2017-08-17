Requests
========

==================
What Is a Request?
==================
A request represents a placeholder for some work to be done.

For example, we could have a request that specifies we want to run a certain
type of job.  Or we could have a request that indicates we want to run a specific flow.

=============
Why Requests?
=============
Requests can help us in several ways:

#. They help us avoid duplicate work by allowing programs to request work,
   without having to know whether that work has already been done. For example:

   #. A computational chemistry flow makes a request, 'request_a', to generate
      conformers for a molecule.

   #. We process request_a and create a conformer job.

   #. We execute the job and save the job's output.
      
   #. We save a reference to the job in request_a and mark request_a as
      completed.
   
   #. Later, a second flow makes a request, 'request_b', to generate conformers
      for the same molecule as before. 

   #. When we process request_b, we see that the job it would have created
      alread exists.  Instead of creating the job again we just save a
      reference to it in request_b and mark request_b as completed.

#. They allow programs to request computations without having to know how those
   computations should be executed. For example, a computational chemistry 
   program can request that we generate conformers. All the chemistry program
   needs to know is how to create a request. We can handle the details of how
   to process that request in a separate system.

#. They allow us to group and organize computations. We can tag requests
   so that multiple projects can share the same database.

================
Request Lifecyle
================

This is the typical request lifecycle.

#. Some program calls a MissionControl command.  For example, a chemistry
   program calls the 
   :class:`mc.houston.subcommands.create_job_requests.Subcommand` command via 
   :doc:`houston` .

# The command creates :class:`mc.db.models.Request` instances in the
  MissionControl DB.

#. A program calls another MissionControl command that updates the request.
   For example, a program calls the
   :class:`mc.houston.subcommands.process_executed_job_dirs.Subcommand` command,
   which ingests executed job dirs and updates requests that refer to the
   job dirs' related jobs.

This is just one lifecycle. Some users may want to use the Request independently
of other parts of MissionControl. For example, to track whether an order has
already been sent to a chemistry robot.

=================
The Request Model
=================
Refer to the documentation for the :class:`mc.db.models.Request` class.


=================
Querying Requests
=================
Usually you won't have to query for requests directly. But if you do, you may
find it useful to look at the code and documentation for the
:class:`mc.utils.selectors.basic_request_selector.BasicRequestSelector` class.
