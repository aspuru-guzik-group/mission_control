Installation
============

=================
Option 1: Use Pip
=================
This is the easiest way to install to MissionControl. This is probably what you want.

#. Run:
   ::

     pip install mission_control

=============================
Option 2: Install from Source
=============================

#. Clone the MissionControl repo.

#. Navigate to the cloned repo and run:
   ::

     python setup.py develop


========================
Verify Your Installation
========================

#. Verify your installation by running:
   ::

     python -m mc.houston.cli sanity_check
