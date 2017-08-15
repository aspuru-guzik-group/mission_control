.. MissionControl documentation master file, created by
   sphinx-quickstart on Tue Jun 27 15:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: Introduction to MissionControl (workflow software)

==============
MissionControl
==============
A library for creating jobs and workflows.

=========================
Is MissionControl for me?
=========================

Features include:

#. A job creation framework: create small job modules that define how to build job directories.
#. A task sequencing framework: create flows that consist of small tasks.
#.  A data storage framework: store data from parsed jobs, and query it to provide inputs to other jobs.
#. A SqlAlchemy backend lets you run MissionControl with a variety of SQL environments. You can use file-based sqlite databases or server-based databases like Postgresql or MySql.
#.  Dynamic flows: you can create flows that modify themselves or create new flows based on what happens during execution.
#. Modular: you can customize how MissionControl works.

========
Overview
========
MissionControl is a collection of several components.

Jobs
====
You can build parameterized jobs that you can run locally or on a cluster.

For example, you can build a set of jobs that run a computational chemistry 
model with different parameters.

See `user_guide/jobs`.

Work Flows
=====
You can use define and run task-based workflows.

For example, you can define a workflow that first runs a model, and then runs
an analysis job on the model outputs.

See `user_guide/flows`.

Requests
========
You can track whether a job or flow has already been
executed with a given set of inputs. This prevents the running of duplicate
computations.

See `user_guide/requests`.

EntityDB
========
You can store data in a generalized
`Entity-Attribute-Value Database
<https://en.wikipedia.org/wiki/Entity%E2%80%93attribute%E2%80%93value_model>`_.

This database can store results from parsed jobs that you
want to use as inputs for future jobs.

See `user_guide/entity_db`.


Combining Components
====================
You can use these components in any combination you choose. Some people use 
MissionControl only to build computing cluster jobs. Others only use
MissionControl to define and run flows.

MissionControl allows you to use components on their own, or in concert.

===============
Getting Started
===============

.. toctree::
   :maxdepth: 1

   quickstart
   user_guide
   examples


=================
API Documentation
=================

.. toctree::
   :maxdepth: 1

   api

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
