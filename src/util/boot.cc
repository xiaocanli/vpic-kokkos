#include "util.h"
#include "git_version.h"

#include <iostream>
#include <omp.h>
#include "../vpic/kokkos_tuning.hpp"

double _boot_timestamp = 0;

double
uptime( void ) {
  double local_time = wallclock(), time_sum;
  mp_allsum_d( &local_time, &time_sum, 1 );
  return time_sum/(double)world_size - _boot_timestamp;
}

void
boot_services( int * pargc,
               char *** pargv )
{
    // Start up the checkpointing service.  This should be first.

    boot_checkpt( pargc, pargv );

    // Start up the threads.  Note that some MPIs will bind threads to
    // cores if threads are booted _after_ MPI is initialized.  So we
    // start up the pipeline dispatchers _before_ starting up MPI.

    // FIXME: The thread utilities should take responsibility for
    // thread-core affinity instead of leaving this to chance.

    serial.boot( pargc, pargv );
    thread.boot( pargc, pargv );

    // Boot up the communications layer
    // See note above about thread-core-affinity

    boot_mp( pargc, pargv );

    // Set the boot_timestamp
    mp_barrier();
    _boot_timestamp = 0;
    _boot_timestamp = uptime();

    // Do some argument clean up to try and save us from ourselves
    if (world_rank == 0)
    {
        // If we have arguments VPIC doesn't understand, notify the user and pass
        // them to Kokkos
        if (*pargc >= 2)
        {
            std::cerr << "Passing the following non-standard args to kokkos:" << std::endl;
            for (int i = 1; i < *pargc; i++)
            {
                std::cerr << i << ") " <<  (*pargv)[i] << std::endl;
            }
            std::cerr << std::endl;
        }

        int nThreads = 0;
#pragma omp parallel
        {
            nThreads = omp_get_num_threads();
        }

        if (nThreads != thread.n_pipeline)
        {
            std::cerr << "omp_get_num_threads != n_pipeline ";
            std::cerr << "(" << nThreads << " != " << thread.n_pipeline << ")" << std::endl;
        }

        std::cerr << "-> Setting omp_get_num_threads to be tpp! " << thread.n_pipeline << std::endl;
        std::cerr << std::endl;
    }

    omp_set_num_threads(thread.n_pipeline);

    // TODO: move this into a specific run outputter class
    if (world_rank == 0)
    {
        std::cout << "######### Build Details ##########" << std::endl;

        std::cout << "# VPIC Git Hash: "  << GIT_REVISION << std::endl;
        std::cout << "# Built on: "  << BUILD_TIMESTAMP << std::endl;
        std::cout << "## CMake options ##" << std::endl;
# if defined ENABLE_INTEGRATED_TESTS
        std::cout << "# ENABLE_INTEGRATED_TESTS: ON" << std::endl;
# else
        std::cout << "# ENABLE_INTEGRATED_TESTS: OFF" << std::endl;
#endif
# if defined ENABLE_UNIT_TESTS
        std::cout << "# ENABLE_UNIT_TESTS: ON" << std::endl;
# else
        std::cout << "# ENABLE_UNIT_TESTS: OFF" << std::endl;
#endif
# if defined USE_V4_ALTIVEC
        std::cout << "# USE_V4_ALTIVEC: ON" << std::endl;
# else
        std::cout << "# USE_V4_ALTIVEC: OFF" << std::endl;
#endif
# if defined USE_V4_PORTABLE
        std::cout << "# USE_V4_PORTABLE: ON" << std::endl;
# else
        std::cout << "# USE_V4_PORTABLE: OFF" << std::endl;
#endif
# if defined USE_V4_SSE
        std::cout << "# USE_V4_SSE: ON" << std::endl;
# else
        std::cout << "# USE_V4_SSE: OFF" << std::endl;
#endif
# if defined ENABLE_OPENSSL
        std::cout << "# ENABLE_OPENSSL: ON" << std::endl;
# else
        std::cout << "# ENABLE_OPENSSL: OFF" << std::endl;
#endif
# if defined ENABLE_KOKKOS_OPENMP
        std::cout << "# ENABLE_KOKKOS_OPENMP: ON" << std::endl;
# else
        std::cout << "# ENABLE_KOKKOS_OPENMP: OFF" << std::endl;
#endif
# if defined ENABLE_KOKKOS_CUDA
        std::cout << "# ENABLE_KOKKOS_CUDA: ON" << std::endl;
# else
        std::cout << "# ENABLE_KOKKOS_CUDA: OFF" << std::endl;
#endif
# if defined BUILD_INTERNAL_KOKKOS
        std::cout << "# BUILD_INTERNAL_KOKKOS: ON" << std::endl;
# else
        std::cout << "# BUILD_INTERNAL_KOKKOS: OFF" << std::endl;
#endif
# if defined VPIC_DUMP_ENERGIES
        std::cout << "# VPIC_DUMP_ENERGIES: ON" << std::endl;
# else
        std::cout << "# VPIC_DUMP_ENERGIES: OFF" << std::endl;
#endif
# if defined VPIC_ENABLE_AUTO_TUNING
        std::cout << "# VPIC_ENABLE_AUTO_TUNING: ON" << std::endl;
# else
        std::cout << "# VPIC_ENABLE_AUTO_TUNING: OFF" << std::endl;
#endif
        std::cout << "# NOTE: Auto tuning may have overwritten what CMake set these as." << std::endl;
# if defined VPIC_ENABLE_HIERARCHICAL
        std::cout << "# VPIC_ENABLE_HIERARCHICAL: ON" << std::endl;
# else
        std::cout << "# VPIC_ENABLE_HIERARCHICAL: OFF" << std::endl;
#endif
# if defined VPIC_ENABLE_TEAM_REDUCTION
        std::cout << "# VPIC_ENABLE_TEAM_REDUCTION: ON" << std::endl;
# else
        std::cout << "# VPIC_ENABLE_TEAM_REDUCTION: OFF" << std::endl;
#endif
# if defined VPIC_ENABLE_VECTORIZATION
        std::cout << "# VPIC_ENABLE_VECTORIZATION: ON" << std::endl;
# else
        std::cout << "# VPIC_ENABLE_VECTORIZATION: OFF" << std::endl;
#endif
# if defined VPIC_ENABLE_ACCUMULATOR
        std::cout << "# VPIC_ENABLE_ACCUMULATOR: ON" << std::endl;
# else
        std::cout << "# VPIC_ENABLE_ACCUMULATOR: OFF" << std::endl;
#endif
        std::cout << "# Using sort method: " << EXPAND_AND_STRINGIFY(SORT) << std::endl;

        std::cout << "######### End Build Details ########" << std::endl;
        std::cout << std::endl; // blank line
    }

    Kokkos::initialize( *pargc, *pargv );
}

// This operates in reverse order from boot_services

void
halt_services( void ) {
  _boot_timestamp = 0;
  halt_mp();
  thread.halt();
  serial.halt();
  halt_checkpt();
}

