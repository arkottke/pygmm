.. highlight:: shell

============
Installation
============

Prior to using `pygmm`, Python and the following dependencies need to be
installed:

* matplotlib -- used for plotting

* numpy -- fast vector operations

`pygmm` supports both Python 2.7 and Python 3.

Linux
-----

Install `pygmm` dependencies is best accomplished with a package manager. On Arch
Linux this can be accomplished with::

    pacman -S python-numpy python-matplotlib pip

Windows and OS-X
----------------

On Windows, installing matplotlib and numpy can be simplified by using
Miniconda3. Miniconda3 has installers for `Windows 32-bit`_, `Windows 64-bit`_,
and `OS-X`_.

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


Installing `pygmm`
------------------

After the dependencies have been installed, install or upgrade `pygmm` using pip::

  pip install --upgrade pygmm

This command can be re-run later to upgrade `pygmm` to the latest version.
