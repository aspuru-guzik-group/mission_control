Components
==========

MissionControl is a collection of several components.

====
Jobs
====

MissionControl's job utilities allow you to build parameterized jobs. You can
run these jobs locally, or on a computing cluster.

For example, you can build a set of jobs that run a computational chemistry 
model with different parameters.

See :doc:`jobs` .

=========
WorkFlows
=========
You can use define and run task-based workflows.

For example, you can define a workflow that runs a model and then analyzes the
outputs.

See :doc:`flows` and :doc:`tasks` .

========
Requests
========

You can track whether a job or flow has already been
executed with a given set of inputs. This helps you avoid running duplicate
computations.

See :doc:`requests` .

========
EntityDB
========
You can store data in a generalized
`Entity-Attribute-Value Database
<https://en.wikipedia.org/wiki/Entity%E2%80%93attribute%E2%80%93value_model>`_.

This database can store results from parsed jobs that you
want to use as inputs for future jobs.

See :doc:`entity_db` .

=======
Houston
=======
Houston is a utility for coordinating and configuring MissionControl components.

It makes it easier to write programs that interact with several parts of
MissionControl.

See :doc:`houston` .


====================
Combining Components
====================
You can use these components in any combination you choose. Some people use 
MissionControl only to build computing cluster jobs. Others only use
MissionControl to define and run flows.

MissionControl allows you to use components on their own, or in concert.
