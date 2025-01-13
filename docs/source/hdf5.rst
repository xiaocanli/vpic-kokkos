HDF5 output for VPIC-kokkos
===========================

How to compile and run?
-----------------------

To compile and run a deck, we need to compile ``VPIC-kokkos`` with HDF5
enabled first. Before the merge to
``https://github.com/lanl/vpic-kokkos``, you need to get it from a fork.
Here, I use Perlmutter@NERSC machine as an example.

.. code:: sh

   git clone https://github.com/xiaocanli/vpic-kokkos
   cd vpic-kokkos
   git checkout xiaocanli/hdf5_output
   mkdir build
   cd build
   ../arch/Perlmutter_GPU

How to change the deck to use HDF5?
-----------------------------------

- For any deck, add these to near the top of the deck \```c++ // Whether
  to use HDF5 format for dumping fields and hydro #define DUMP_WITH_HDF5

  #ifdef DUMP_WITH_HDF5 // Deck only works if VPIC was build with HDF
  support. Check for that: #ifndef VPIC_ENABLE_HDF5 #error
  “VPIC_ENABLE_HDF5” is required #endif #endif \``\`

- Add these to the end of ``begin_initialization { }`` \```c++ #ifdef
  DUMP_WITH_HDF5 // For writing XDMF file when using HDF5 dump
  field_interval = global->fields_interval; hydro_interval =
  global->ehydro_interval;

  ::

     enable_hdf5_dump();

  #endif
  :literal:`Remember to load \`HDF5\` module before compiling and running. On Perlmutter,`\ sh
  module load cray-hdf5-parallel \``\` Then, we can compile and run the
  deck as usual.

How to analyze the data?
------------------------

The code will output fields and hydro data to ``field_hdf5`` and
``hydro_hdf5``, respectively. Note that electron hydro data and ion
hydro data are saved in different files. You can check what is included
in the ``.h5`` files using ``h5dump`` after loading the HDF5 modules. To
learn more about these command-line tools, please check out `HDF5
TUTORIAL: COMMAND-LINE TOOLS FOR VIEWING HDF5
FILES <https://support.hdfgroup.org/HDF5/Tutor/cmdtoolview.html>`__. The
whole `HDF5 TUTORIAL <https://support.hdfgroup.org/HDF5/Tutor/>`__ is
very helpful if you are not familiar with HDF5. For example,
``h5dump -H fields_0.h5`` gives

::

   HDF5 "fields_0.h5" {
   GROUP "/" {
      GROUP "Timestep_0" {
         DATASET "cbx" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
         DATASET "cby" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
         DATASET "cbz" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
         DATASET "ex" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
         DATASET "ey" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
         DATASET "ez" {
            DATATYPE  H5T_IEEE_F32LE
            DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
         }
      }
   }
   }

``cbx``, ``cby``, ``cbz`` are magnetic fields and ``ex``, ``ey``, ``ez``
are electric fields. Similarly, the hydro files include ``jx``, ``jy``,
and ``jz`` for current density, ``ke`` for kinetic energy density,
``px``, ``py``, and ``pz`` for momentum density, ``rho`` for charge
density, ``txx``, ``txy``, ``tyy``, ``tyz``, ``tzx``, and ``tzz`` for
stress tensor. These are similar to the data before ``translate``. -
``jx/rho`` will be ``vex`` or ``vix`` - ``|rho|`` will be ``ne`` or
``ni`` - ``px/|rho|/particle_mass`` with be the commonly used ``uex`` or
``uix``. ``particle_mass`` is 1 for electrons or ``mi_me`` for ions. -
``txx - (jx/rho)*px`` will be ``pexx`` or ``pixx``.
``txy - (jx/rho)*py`` will be ``pexy`` or ``pixy``.
``txy - (jy/rho)*px`` will be ``peyx`` or ``piyx``. Similar for other
pressure tensor components.

You can use
`quick_check_vpic <https://github.com/xiaocanli/quick_check_vpic>`__ to
check the data. For detailed analysis using the HDF5 files, you need to
write your own scripts. For example,

.. code:: py

   import h5py
   import matplotlib.pylab as plt
   import numpy as np

   field_interval = 200
   hydro_interval = 200
   tframe = 1

   # read the fields
   tindex = tframe * field_interval
   fname = "fields_hdf5/T." + str(tindex) + "/fields_" + str(tindex) + ".h5"
   with h5py.File(fname, "r") as fh:
       group = fh["Timestep_" + str(tindex)]
       # list all the fields
       print(group.keys())
       # read Bx for example
       Bx = group["cbx"][:, :, :]  # read all the data
       print(Bx.shape)
       plt.imshow(Bx[:, 0, :].T, origin="lower")
       plt.show()
       # alternatively, we can read a subset conviniently
       dset = group["cbx"]
       Bx = dset[480:, 0, :320]  # read a subset
       plt.imshow(Bx.T, origin="lower")
       plt.show()

It reads and plot :math:`B_x`. Reading other variables is similar. After
reading into the memory, the rest will be the same. You can still use
most of your analysis code. If the data is really large (e.g., in 3D
simulations), we use ParaView or VisIt to visualize the data. In
``fields_hdf5/`` and ``hydro_hdf5/``, there are files ending with
``.xdmf``, which can be loaded into ParaView or VisIt directly.
