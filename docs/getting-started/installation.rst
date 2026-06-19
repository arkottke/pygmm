============
Installation
============

pyGMM supports Python 3.10 and later versions. We recommend using the latest stable version of Python.

.. grid:: 1 1 2 2
    :gutter: 3

    .. grid-item-card:: 🚀 Quick Install
        :class-title: card-title-center

        For most users, install with pip:

        .. code-block:: bash

           pip install pygmm

    .. grid-item-card:: 🔧 Development Install
        :class-title: card-title-center

        For contributors and developers:

        .. code-block:: bash

           git clone https://github.com/arkottke/pygmm.git
           cd pygmm
           pip install -e .

Requirements
============

pyGMM requires the following packages:

.. grid:: 1 2 3 3
    :gutter: 2

    .. grid-item::
        :columns: 12 6 4 4

        **Core Dependencies**

        - Python 3.10+
        - NumPy
        - SciPy
        - Matplotlib

    .. grid-item::
        :columns: 12 6 4 4

        **Optional Dependencies**

        - pandas (for data analysis)
        - Jupyter (for notebooks)

    .. grid-item::
        :columns: 12 12 4 4

        **Development Dependencies**

        - pytest (testing)
        - sphinx (documentation)
        - ruff (linting)

Installation Methods
====================

.. tab-set::

    .. tab-item:: pip (Recommended)

        Install the latest stable release from PyPI:

        .. code-block:: bash

           pip install pygmm

        To upgrade to the latest version:

        .. code-block:: bash

           pip install --upgrade pygmm

    .. tab-item:: conda

        Install from conda-forge:

        .. code-block:: bash

           conda install -c conda-forge pygmm

        Or using mamba:

        .. code-block:: bash

           mamba install -c conda-forge pygmm

    .. tab-item:: Development

        For the latest development version:

        .. code-block:: bash

           pip install git+https://github.com/arkottke/pygmm.git

        Or clone the repository:

        .. code-block:: bash

           git clone https://github.com/arkottke/pygmm.git
           cd pygmm
           pip install -e .[test,docs]

Virtual Environments
=====================

We strongly recommend using virtual environments to avoid conflicts:

.. tab-set::

    .. tab-item:: venv

        .. code-block:: bash

           python -m venv pygmm-env
           source pygmm-env/bin/activate  # On Windows: pygmm-env\Scripts\activate
           pip install pygmm

    .. tab-item:: conda

        .. code-block:: bash

           conda create -n pygmm-env python=3.11
           conda activate pygmm-env
           pip install pygmm

    .. tab-item:: uv (Fast)

        Using the modern `uv` package manager:

        .. code-block:: bash

           uv venv pygmm-env
           source pygmm-env/bin/activate  # On Windows: pygmm-env\Scripts\activate
           uv pip install pygmm

Verification
============

Test your installation:

.. code-block:: python

   import pygmm
   print(f"pyGMM version: {pygmm.__version__}")

   # Quick test
   scenario = pygmm.Scenario(mag=6.0, dist_rup=10, v_s30=760)
   model = pygmm.CampbellBozorgnia2014()
   result = model(scenario)
   print("Installation successful!")

Troubleshooting
===============

.. dropdown:: Common Issues
   :class-title: sd-bg-info sd-text-white

   **Import Error**

   If you get import errors, ensure you've activated the correct environment:

   .. code-block:: bash

      which python  # Should point to your virtual environment
      pip list      # Check if pygmm is installed

   **Missing Dependencies**

   Install missing dependencies:

   .. code-block:: bash

      pip install numpy scipy matplotlib

   **Permission Errors**

   Use ``--user`` flag or virtual environments:

   .. code-block:: bash

      pip install --user pygmm

.. dropdown:: Platform-Specific Notes
   :class-title: sd-bg-secondary sd-text-white

   **Windows**

   - Use Command Prompt or PowerShell
   - Consider Windows Subsystem for Linux (WSL) for better compatibility

   **macOS**

   - Install Xcode command line tools: ``xcode-select --install``
   - Consider using Homebrew for Python: ``brew install python``

   **Linux**

   - Most distributions include Python 3.10+
   - Install development headers if building from source

Getting Help
============

If you encounter issues:

1. Check the :doc:`examples/index` for common usage patterns
2. Search existing `GitHub Issues <https://github.com/arkottke/pygmm/issues>`_
3. Create a new issue with detailed information about your problem
