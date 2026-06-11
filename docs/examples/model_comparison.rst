================
Model Comparison
================

This example shows how to compare predictions from multiple ground motion models.

Comparing Multiple Models
=========================

pyGMM makes it easy to compare different ground motion prediction equations:

.. code-block:: python

   import pygmm
   import numpy as np
   import matplotlib.pyplot as plt

   # Define scenario
   scenario = pygmm.Scenario(
       mag=6.5,
       dist_rup=20.0,
       v_s30=760,
       mechanism='strike_slip'
   )

   # Initialize multiple models
   models = {
       'CB14': pygmm.CampbellBozorgnia2014(),
       'ASK14': pygmm.AbrahamsonSilvaKamai2014(),
       'BSSA14': pygmm.BooreStewartSeyhanAtkinson2014(),
       'CY14': pygmm.ChiouYoungs2014(),
   }

   # Calculate predictions for each model
   results = {}
   for name, model in models.items():
       ln_sa, ln_std = model(scenario)
       results[name] = {
           'sa': np.exp(ln_sa),
           'periods': model.periods
       }

Plotting Comparison
===================

Create a comparison plot:

.. code-block:: python

   plt.figure(figsize=(12, 8))

   colors = ['blue', 'red', 'green', 'orange']

   for i, (name, data) in enumerate(results.items()):
       plt.loglog(data['periods'], data['sa'],
                  color=colors[i], label=name, linewidth=2)

   plt.xlabel('Period (s)')
   plt.ylabel('Spectral Acceleration (g)')
   plt.title('Ground Motion Model Comparison')
   plt.legend()
   plt.grid(True, alpha=0.3)
   plt.xlim(0.01, 10)
   plt.show()

Statistical Analysis
====================

Analyze the variability between models:

.. code-block:: python

   # Extract SA values at common periods
   common_periods = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

   model_data = []
   for period in common_periods:
       period_data = []
       for name, data in results.items():
           # Find closest period
           idx = np.argmin(np.abs(data['periods'] - period))
           period_data.append(data['sa'][idx])
       model_data.append(period_data)

   # Calculate statistics
   means = [np.mean(data) for data in model_data]
   stds = [np.std(data) for data in model_data]
   coeffs_of_var = [std/mean for mean, std in zip(means, stds)]

   # Display results
   print("Period (s) | Mean SA (g) | Std Dev | CoV")
   print("-" * 40)
   for i, period in enumerate(common_periods):
       print(f"{period:8.2f} | {means[i]:10.4f} | {stds[i]:7.4f} | {coeffs_of_var[i]:5.3f}")

.. tip::

   Different models may have different strengths depending on the scenario.
   Consider the applicable ranges and uncertainties when selecting models.
