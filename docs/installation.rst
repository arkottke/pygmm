.. highlight:: shell

============
Installation
============

Prior to using pyGMM, Python and a number of packages need to be installed. In
addition to Python, the following packages need to be installed:

* matplotlib -- used for plotting

* numpy -- fast vector operations

Install Python dependencies is best accomplished with a package manager. On
Windows or OS-X, I recommend using Miniconda3. On Linux, the package manager
is preferred.

Miniconda3 has installers for `Windows 32-bit`_, `Windows 64-bit`_, and `OS-X`_.
While pyGMM has support Python27, Python3 is preferred.

.. _Windows 32-bit: http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe
.. _Windows 64-bit: http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe
.. _OS-X: http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

After the installer is finished, install the required dependencies by opening a
terminal. On Windows, this is best accomplished with ``Windows Key + r``, enter
``cmd``. Next enter the following command::

  conda install --yes pip setuptools numpy matplotlib

On Windows, the text can copied and pasted if *Quick Edit* mode is enabled. To
enable this feature, right click on the icon in the upper left portion of the
window, and select *Properties*, and then check the *Quick Edit Mode* check box
within the *Edit Options* group. Copy the text, and then paste it by click the
right mouse button.

Now that the dependencies have been installed, install or upgrade pyGMM using
pip::

  pip install --upgrade pygmm

This command can be re-run later to upgrade pyGMM to the latest version.
