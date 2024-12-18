.. SWIRL documentation master file, created by
   sphinx-quickstart on Sun Nov 12 13:23:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/SWIRL_logo.png
   :width: 300
   
SWIRL - *S*\tructural *W*\ind-borne debris *I*\mpact *R*\isk assessment *L*\ibrary
==================================================================================

**Author**: `Brian Doran Giffin <https://github.com/bdgiffin>`_ (`brian.giffin@okstate.edu <mailto:brian.giffin@okstate.edu>`_)

Wind-borne debris is a significant contributor to structural damage during high-intensity wind events, but existing methods for estimating debris impact loads on structures are relatively limited. These limitations stem from inherent uncertainties and lack of knowledge regarding the characterization of combined wind and debris loads, as well as a lack of computational modeling strategies for representing wind-borne debris impacts and their effects on structures.

The **S**\tructural **W**\ind-borne debris **I**\mpact **R**\isk assessment **L**\ibrary (**SWIRL**) is a physics-based fluid-structure-debris modeling framework intended to be coupled with the open source structural analysis software `OpenSEES <https://opensees.berkeley.edu>`_. SWIRL was originally developed to investigate and quantify of the extent to which wind-borne debris impact contributes to structural damage and collapse. Within this modeling framework, flying debris is represented through discrete realizations of debris trajectories and impacts, with the nonlinear transient dynamic behavior of the structure of interest modeled using OpenSees. Collisions between debris and the structure are resolved through a penalty-based contact enforcement strategy. Parametric vortex models are used to represent the wind field and to determine wind pressures acting on both the structure and the debris.

The core functionality of SWIRL is written in C++, with a fully supported Python API. SWIRL is compatible with OpenSeesPy, and can be used in Python workflows that can be run using `quoFEM <https://simcenter.designsafe-ci.org/research-tools/quofem-application/>`_. Several simple :ref:`Examples` are provided to demonstrate how SWIRL may be utilized to conduct a variety of risk assessment analyses.

Information regarding the :ref:`Installation` of SWIRL and additional :ref:`Resources` may be found within the accompanying documentation pages.

.. note::

   The documentation for SWIRL is under active development. Please check back periodically for updates.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home <self>
   installation
   examples
   resources

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
