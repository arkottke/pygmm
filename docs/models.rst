.. _models:

Ground motion models
====================

This part of the documentation covers all the interfaces of pyGMM.

Generic Interface
-----------------

The interfaces for models have been simplified to use same parameter names and
values where possible.


+--------------+--------------+
| Abbreviation | Name         |
+==============+==============+
| u            | Unspecified  |
+--------------+--------------+
| ss           | Strike-slip  |
+--------------+--------------+
| ns           | Normal slip  |
+--------------+--------------+
| rs           | Reverse slip |
+--------------+--------------+


All model share the following interface:

.. autoclass:: pygmm.model.Model
    :members:



Specific Models
---------------

.. automodule:: pygmm
    :members:
    :special-members:
