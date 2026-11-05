Scenarios
=========

A :class:`~pygmm.model.Scenario` describes the earthquake and the site. It is a
dictionary with a fixed set of keys, so a misspelled parameter raises an error
rather than being ignored.

.. code-block:: python

   import pygmm

   scenario = pygmm.Scenario(
       mag=6.5,
       dist_rup=20.0,
       dist_jb=20.0,
       dist_x=20.0,
       dip=90.0,
       mechanism="SS",
       v_s30=760.0,
   )

Parameters
----------

The keys, with their units, are listed in the documentation of
:class:`~pygmm.model.Scenario`. Distances are in km, depths in km, angles in
degrees, :math:`V_{s30}` in m/s and the magnitude is moment magnitude.

Each model uses a subset of the keys. The parameters of a model are listed in its
``PARAMS`` attribute, and a model raises ``ValueError`` for a required parameter
that is missing, for example ``('dip', 'is a required parameter')``. Keys that a
model does not use are ignored.

Fault mechanism
---------------

.. list-table::
   :header-rows: 1

   * - Abbreviation
     - Name
   * - ``U``
     - Unspecified
   * - ``SS``
     - Strike-slip
   * - ``NS``
     - Normal slip
   * - ``RS``
     - Reverse slip

Each model defines the mechanism for itself, and some accept a different set of
values.

Copies
------

``Scenario.copy_with(**kwds)`` returns a scenario with some keys changed, which
is how a parameter is varied across models or across a range of values; see
:doc:`../examples/example-01`.
