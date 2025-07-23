=======================
Probabilistic Analysis
=======================

This example shows how to perform probabilistic seismic hazard analysis using pyGMM.

Monte Carlo Simulation
======================

Propagate uncertainty through ground motion calculations:

.. code-block:: python

   import pygmm
   import numpy as np
   import matplotlib.pyplot as plt
   from scipy import stats

   # Define base scenario
   base_scenario = pygmm.Scenario(
       mag=6.5,
       dist_rup=20,
       v_s30=760,
       mechanism='strike_slip'
   )

   # Define uncertainty in input parameters
   n_samples = 1000
   np.random.seed(42)  # For reproducibility

   # Sample magnitude with uncertainty
   mag_samples = np.random.normal(6.5, 0.2, n_samples)

   # Sample distance with uncertainty
   dist_samples = np.random.lognormal(np.log(20), 0.3, n_samples)

   # Sample Vs30 with uncertainty
   vs30_samples = np.random.lognormal(np.log(760), 0.4, n_samples)

Running Monte Carlo
====================

Perform the analysis:

.. code-block:: python

   model = pygmm.CampbellBozorgnia2014()

   # Store results
   pga_results = []
   sa_1s_results = []

   for i in range(n_samples):
       # Create scenario for this realization
       scenario = pygmm.Scenario(
           mag=mag_samples[i],
           dist_rup=dist_samples[i],
           v_s30=vs30_samples[i],
           mechanism='strike_slip'
       )

       # Calculate ground motion
       ln_sa, ln_std = model(scenario)
       sa = np.exp(ln_sa)

       # Add aleatory uncertainty
       epsilon = np.random.normal(0, 1, len(ln_sa))
       sa_with_aleatory = sa * np.exp(epsilon * np.exp(ln_std))

       # Store results
       pga_results.append(sa_with_aleatory[model.INDEX_PGA])

       # SA at 1.0 second
       idx_1s = np.argmin(np.abs(model.periods - 1.0))
       sa_1s_results.append(sa_with_aleatory[idx_1s])

   pga_results = np.array(pga_results)
   sa_1s_results = np.array(sa_1s_results)

Statistical Analysis
====================

Analyze the results:

.. code-block:: python

   # Calculate statistics
   print("Probabilistic Ground Motion Results")
   print("=" * 40)
   print(f"PGA Statistics:")
   print(f"  Mean: {np.mean(pga_results):.3f} g")
   print(f"  Median: {np.median(pga_results):.3f} g")
   print(f"  Standard Deviation: {np.std(pga_results):.3f} g")
   print(f"  84th Percentile: {np.percentile(pga_results, 84):.3f} g")
   print(f"  16th Percentile: {np.percentile(pga_results, 16):.3f} g")

   print(f"\nSA(1.0s) Statistics:")
   print(f"  Mean: {np.mean(sa_1s_results):.3f} g")
   print(f"  Median: {np.median(sa_1s_results):.3f} g")
   print(f"  Standard Deviation: {np.std(sa_1s_results):.3f} g")

Visualization
=============

Create probability plots:

.. code-block:: python

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

   # PGA histogram
   ax1.hist(pga_results, bins=50, density=True, alpha=0.7,
            edgecolor='black', label='Simulated')

   # Fit lognormal distribution
   sigma, loc, scale = stats.lognorm.fit(pga_results, floc=0)
   x = np.linspace(0.001, max(pga_results), 1000)
   ax1.plot(x, stats.lognorm.pdf(x, sigma, loc=loc, scale=scale),
            'r-', linewidth=2, label='Lognormal fit')

   ax1.set_xlabel('PGA (g)')
   ax1.set_ylabel('Probability Density')
   ax1.set_title('PGA Distribution')
   ax1.legend()
   ax1.grid(True, alpha=0.3)

   # SA(1.0s) histogram
   ax2.hist(sa_1s_results, bins=50, density=True, alpha=0.7,
            edgecolor='black', label='Simulated')

   sigma, loc, scale = stats.lognorm.fit(sa_1s_results, floc=0)
   x = np.linspace(0.001, max(sa_1s_results), 1000)
   ax2.plot(x, stats.lognorm.pdf(x, sigma, loc=loc, scale=scale),
            'r-', linewidth=2, label='Lognormal fit')

   ax2.set_xlabel('SA(1.0s) (g)')
   ax2.set_ylabel('Probability Density')
   ax2.set_title('SA(1.0s) Distribution')
   ax2.legend()
   ax2.grid(True, alpha=0.3)

   plt.tight_layout()
   plt.show()

Hazard Curves
=============

Calculate exceedance probabilities:

.. code-block:: python

   # Define intensity levels
   pga_levels = np.logspace(-3, 0, 50)  # 0.001 to 1.0 g

   # Calculate exceedance probabilities
   exceedance_probs = []
   for level in pga_levels:
       prob = np.mean(pga_results > level)
       exceedance_probs.append(prob)

   # Plot hazard curve
   plt.figure(figsize=(10, 6))
   plt.loglog(pga_levels, exceedance_probs, 'b-', linewidth=2)
   plt.xlabel('PGA (g)')
   plt.ylabel('Probability of Exceedance')
   plt.title('PGA Hazard Curve')
   plt.grid(True, alpha=0.3)
   plt.show()

   # Find PGA at specific probability levels
   target_probs = [0.1, 0.02, 0.002]  # 10%, 2%, 0.2% probability

   print("\nPGA at specific probability levels:")
   for prob in target_probs:
       idx = np.argmin(np.abs(np.array(exceedance_probs) - prob))
       pga_at_prob = pga_levels[idx]
       print(f"  {prob*100:4.1f}% probability: {pga_at_prob:.3f} g")

.. note::

   This is a simplified example. Real probabilistic seismic hazard analysis
   involves more complex source characterization and integration over all
   possible earthquake scenarios.
