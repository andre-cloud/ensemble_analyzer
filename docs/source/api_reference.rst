API Reference
=============

This section documents the main classes and functions of **Ensemble Analyzer**.

.. module:: ensemble_analyzer

Data Structures
---------------
Fundamental classes defining the molecule object and associated data structures.

.. automodule:: ensemble_analyzer._conformer.conformer
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._conformer.energy_data
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer._conformer.spectral_data
   :members:
   :undoc-members:

Workflow Managers
-----------------
The core engine managing the protocol execution, job orchestration, and pruning workflows.

.. automodule:: ensemble_analyzer._managers.protocol_manager
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._managers.pruning_manager
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._managers.calculation_executor
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer._managers.checkpoint_manager
   :members:
   :undoc-members:

Algorithms
---------------------
Implementation of clustering algorithms (EDM eigenvalues) and thermochemistry (qRRHO).

.. automodule:: ensemble_analyzer.clustering
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._clustering.cluster_manager
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._clustering.cluster_config
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer.rrho
   :members:
   :undoc-members:
   :show-inheritance:

Spectroscopy & Plotting
-----------------------
Modules for vibrational and electronic spectra convolution and graph generation.

.. automodule:: ensemble_analyzer.graph
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer._spectral.comp_vibronic
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._spectral.comp_electronic
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ensemble_analyzer._spectral.graph_default
   :members:
   :undoc-members:

QM Drivers & Parsers
--------------------
Interfaces and parsers for external quantum chemistry software packages.

.. automodule:: ensemble_analyzer._calculators.base
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer._calculators._orca
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer._calculators._gaussian
   :members:
   :undoc-members:

Input/Output & Utilities
------------------------

.. automodule:: ensemble_analyzer.ensemble_io
   :members:
   :undoc-members:

.. automodule:: ensemble_analyzer.io_utils
   :members:
   :undoc-members: