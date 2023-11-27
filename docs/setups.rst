Initial conditions / pre-cooked configurations
==============================================

How to select a pre-cooked compile-time configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When writing your local Makefile simply give the name of the desired SETUP variable as the argument. For example::

   ~/phantom/scripts/writemake.sh disc > Makefile

Alternatively you can override this manually at compile-time::

   make SETUP=disc
   make SETUP=disc setup

Possible values of the SETUP variable are listed below.

Well-maintained setups
~~~~~~~~~~~~~~~~~~~~~~

The most used, well-maintained and supported options for the SETUP variable are:

.. include:: setups-best.rst

All setups
~~~~~~~~~~~
The full list of pre-cooked configurations (taken from `build/Makefile_setups <https://github.com/danieljprice/phantom/blob/master/build/Makefile_setups>`__) is as follows:

.. include:: setups-list.rst

More information
~~~~~~~~~~~~~~~~~
For a fuller explanation of compile-time configuration options in Phantom, see :doc:`Compile-time configuration <config>`.
