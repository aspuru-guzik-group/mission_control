.. MissionControl documentation master file, created by
   sphinx-quickstart on Tue Jun 27 15:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: Introduction to MissionControl (workflow software)

.. image:: _static/mc_logo.png
   :width: 300 px
   :alt: MissionControl (workflow software)

MissionControl is a library for creating jobs and workflows.

====================
Is MissionControl for me?
====================

Some (but not all) of its features include:

* A job creation framework: create small job modules that define how to build job directories.

* A task sequencing framework: create flows that consist of small tasks.

* A data storage framework: store data from parsed jobs, and query it to provide inputs to other jobs.

* A SqlAlchemy backend lets you run MissionControl with a variety of SQL environments. You can use file-based sqlite databases or server-based databases like Postgresql or MySql.

* *Dynamic* flows: you can create flows that modify themselves or create new flows based on what happens during execution.

* Modular: you can customize how MissionControl works.

==============================
High-Level View of MissionControl
==============================
The main components of MissionControl are:

#. Jobs: a framework for creating job directories or workunits that can run on high-performance computing clusters or grid clusters.

#. Requests: a framework for tracking whether various types of requests have been executed for a given set of inputs and parameters.

#. EntityDB: a framework for storing results from jobs and querying these results to provide inputs to other jobs.

#. Flows: a framework for defining and running sequences of tasks and jobs.

You can use these components in any combination you choose. Some people only need the JobModule framework. Others only need the Flow framework.

For people who want to combine components, MissionControl provides utilities to help with this.

Jobs
====

JOB MODULES HERE!

Requests
========

EntityDB
========

Flows
=====


===============
Getting Started
===============

Quickstart
==========

Start with our installation and quickstart tutorials:

.. toctree::
   :maxdepth: 1

   installation
   quickstart

Examples
========

After you complete the quickstart checkout the examples

.. toctree::
   :maxdepth: 1

   examples


====================================
Contributing / Contact / Bug Reports
====================================

Want to see something added or changed? There are many ways to make that a reality! Some ways to get involved are:

* Help us improve the documentation - tell us where you got 'stuck' and improve the install process for everyone.
* Let us know if you need support for a queueing system or certain features.
* Point us to areas of the code that are difficult to understand or use.
* Contribute code! If you are interested in this option, please see our :doc:`contribution guidelines</contributing>`.


===========================
API Documentation
===========================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
