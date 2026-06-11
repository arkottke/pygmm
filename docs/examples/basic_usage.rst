===========
Basic Usage
===========

This example demonstrates the basic usage of pyGMM for calculating ground motion predictions.

Creating a Scenario
====================

The first step is to define an earthquake scenario:

.. code-block:: python

   import pygmm

   # Define earthquake and site parameters
   scenario = pygmm.Scenario(
       mag=7.0,                # Magnitude
       dist_rup=25.0,          # Rupture distance (km)
       v_s30=760,              # Site shear-wave velocity (m/s)
       mechanism='strike_slip', # Fault mechanism
   )

Using a Ground Motion Model
===========================

Next, select and initialize a ground motion prediction equation:

.. code-block:: python

   # Initialize the Campbell & Bozorgnia (2014) model
   gmpe = pygmm.CampbellBozorgnia2014()

   # Calculate ground motion predictions
   ln_sa, ln_std = gmpe(scenario)

   # Convert to linear units (g)
   import numpy as np
   sa = np.exp(ln_sa)
   std = np.exp(ln_std)

Accessing Results
=================

The model returns spectral accelerations at standard periods:

.. code-block:: python

   # Get the periods
   periods = gmpe.periods

   # Peak ground acceleration (PGA)
   pga = sa[gmpe.INDEX_PGA]

   # Spectral acceleration at 1.0 second
   sa_1s_idx = np.argmin(np.abs(periods - 1.0))
   sa_1s = sa[sa_1s_idx]

   print(f"PGA: {pga:.3f} g")
   print(f"SA(1.0s): {sa_1s:.3f} g")

Plotting Results
================

Visualize the response spectrum:

.. code-block:: python

   import matplotlib.pyplot as plt

   plt.figure(figsize=(10, 6))
   plt.loglog(periods, sa)
   plt.xlabel('Period (s)')
   plt.ylabel('Spectral Acceleration (g)')
   plt.title('Response Spectrum')
   plt.grid(True, alpha=0.3)
   plt.show()

.. note::

   This is a basic example. For more complex scenarios involving multiple models
   or uncertainty quantification, see the advanced examples.
