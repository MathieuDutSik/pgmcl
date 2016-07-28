Parallel Geophysical Model Coupling Library
===========================================

This fortran library deals with coupling geophysical models.
It allow efficient interfacing between models with very good speed.
It is much simpler to use than several other library.

Overview
--------

Nowadays, it is relatively common to couple geophysical models (say atmosphere
with wave model) in order to get better forecasts results. This library
provides functionality for efficiently
  * Computing the interpolation arrays needed in parallel
  * Exchanging data between models efficiently with never a global array being created
  * Exchanges can be with blocked MPI calls or unblocked.

References
----------

Several models have been coupled by using this library. See following publications:
  * M. Dutour Sikirić, A. Roland, I. Janeković, I. Tomažić, M. Kuzmić, Coupling of the Regional Ocean Modelling System and Wind Wave Model, Ocean Modelling 72 (2013) 59--73.
  * M. Dutour Sikirić, A. Roland, I. Tomažić, I. Janeković, Hindcasting the Adriatic Sea near-surface motions with a coupled wave-current model, Journal of Geophysical Research - Oceans, 117 (2012) C00J36
  * L. Cavaleri, A. Roland, M. Dutour Sikirić, L. Bertotti, L. Torrisi, On the coupling of COSMO to WAM, Proceedings of ECMWF Workshop on Ocean Waves, 25-27 June 2012, edited by J. Bidlot.

Dependencies
------------

The library uses MPI for the coupling and naturally depends also on the model considered.
