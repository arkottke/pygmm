=============
Site Response
=============

This example demonstrates how to account for site effects in ground motion predictions.

Site Classification
====================

Different models use various site parameters:

.. code-block:: python

   import pygmm
   import numpy as np

   # Define scenarios for different site classes
   scenarios = {
       'Rock (NEHRP A)': pygmm.Scenario(
           mag=7.0, dist_rup=30, v_s30=1500, mechanism='strike_slip'
       ),
       'Firm Soil (NEHRP C)': pygmm.Scenario(
           mag=7.0, dist_rup=30, v_s30=760, mechanism='strike_slip'
       ),
       'Soft Soil (NEHRP D)': pygmm.Scenario(
           mag=7.0, dist_rup=30, v_s30=360, mechanism='strike_slip'
       ),
       'Very Soft Soil (NEHRP E)': pygmm.Scenario(
           mag=7.0, dist_rup=30, v_s30=180, mechanism='strike_slip'
       )
   }

Site Amplification Effects
===========================

Calculate and compare site amplification:

.. code-block:: python

   import matplotlib.pyplot as plt

   # Use Campbell & Bozorgnia (2014) model
   model = pygmm.CampbellBozorgnia2014()

   # Calculate for each site class
   results = {}
   for site_name, scenario in scenarios.items():
       ln_sa, ln_std = model(scenario)
       results[site_name] = {
           'sa': np.exp(ln_sa),
           'periods': model.periods
       }

   # Plot comparison
   plt.figure(figsize=(12, 8))

   colors = ['brown', 'blue', 'orange', 'red']
   for i, (site_name, data) in enumerate(results.items()):
       plt.loglog(data['periods'], data['sa'],
                  color=colors[i], label=site_name, linewidth=2)

   plt.xlabel('Period (s)')
   plt.ylabel('Spectral Acceleration (g)')
   plt.title('Site Effects on Response Spectra')
   plt.legend()
   plt.grid(True, alpha=0.3)
   plt.show()

Site Amplification Factors
===========================

Calculate amplification relative to rock:

.. code-block:: python

   # Use rock as reference
   rock_sa = results['Rock (NEHRP A)']['sa']
   periods = results['Rock (NEHRP A)']['periods']

   # Calculate amplification factors
   print("Site Amplification Factors (relative to rock)")
   print("=" * 50)

   for site_name, data in results.items():
       if site_name == 'Rock (NEHRP A)':
           continue

       amplification = data['sa'] / rock_sa

       # Show amplification at key periods
       key_periods = [0.1, 0.2, 0.5, 1.0, 2.0]
       print(f"\n{site_name}:")
       for period in key_periods:
           idx = np.argmin(np.abs(periods - period))
           print(f"  T = {period}s: {amplification[idx]:.2f}")

Advanced Site Parameters
========================

Some models support additional site parameters:

.. code-block:: python

   # Example with basin depth parameters
   scenario_with_basin = pygmm.Scenario(
       mag=6.5,
       dist_rup=25,
       v_s30=360,
       depth_1_0=0.5,    # Depth to 1.0 km/s (km)
       depth_2_5=2.0,    # Depth to 2.5 km/s (km)
       mechanism='strike_slip'
   )

   # Models that support basin effects
   ask14 = pygmm.AbrahamsonSilvaKamai2014()
   ln_sa, ln_std = ask14(scenario_with_basin)

   print(f"With basin effects included")
   print(f"PGA: {np.exp(ln_sa[ask14.INDEX_PGA]):.3f} g")

.. warning::

   Not all models support all site parameters. Check the model documentation
   for required and optional parameters.
