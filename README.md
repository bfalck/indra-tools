# indra-tools
This is the code repository for the Indra simulations, hosted on [SciServer](http://www.sciserver.org).

Indra is a suite of large-volume cosmological *N*-body simulations. Each of the 384 simulations is computed with the same WMAP7 cosmological parameters and different initial phases, providing excellent statistics of the large-scale features of the distribution of dark matter in the universe. The independent volumes each have 1024<sup>3</sup> dark matter particles in a box of length 1 Gpc/*h*.

The Indra data volumes contain, for each simulation:
- 64 snapshots of particle positions and velocities
- 64 snapshots of FOF and SUBFIND halo catalogs
- 505 time-steps of coarse-gridded Fourier-space density fields

The Indra relational database contains:
- Halo catalog tables for every simulation and snapshot
- Spatial3D library to allow efficient selection of halos and particle data within 3-dimensional shapes

## Getting Started
To access Indra, create an account on [SciServer](http://www.sciserver.org) and join the *Cosmological Simulations* Science Domain. Then, from the SciServer *Dashboard*, navigate to *Compute* and create a container: select the *Cosmological Simulations* compute image and mount all three Indra *Data Volumes*. 

The *Cosmological Simulations* compute image comes with this library pre-installed plus all software on the default *SciServer Essentials* image. 

To install by hand, go to a terminal in a SciServer container (that has the Indra data volumes mounted) and execute:

`pip install git+https://github.com/bfalck/indra-tools.git`



Example notebooks are provided that demonstrate how to:
- use the reading functions: `read_examples.ipynb` (**Start here**)
- query the halo catalog database tables: `database_examples.ipynb`
- compute density fields: `density_field_examples.ipynb`
- selectively (and efficiently) read particles in 3-dimensional shapes such as a sphere, box, or cone: `Shape3D_examples.ipynb`

These notebooks can be found in the main *Indra* data volume from within a compute container. Copy them to your *persistent* home directory to use and edit. 


## Reference
indra-tools is licensed under a 3-clause BSD style license.

We ask that scientific publications that make use of Indra cite the Falck, et al. 2021 [data release paper](https://arxiv.org/abs/2101.03631) (Falck, et al. 2021, MNRAS, 506, 2659).
