include(CheckFortranSourceCompiles)

# create interface target with compiler flags
if(NOT TARGET compile-options)
    add_library(compile-options INTERFACE)
endif()


function(allow_mismatched_mpi_args)
    if(NOT TARGET MPI::MPI_Fortran)
        find_package(MPI REQUIRED)
    endif()

    # save state, if any
    if(DEFINED CMAKE_REQUIRED_LIBRARIES)
        set(TMP_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    endif()
    if(DEFINED CMAKE_REQUIRED_FLAGS)
        set(TMP_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    endif()

    set(CMAKE_REQUIRED_LIBRARIES "MPI::MPI_Fortran")
    string(CONCAT mpi_mismatched_args_prog
        "program main
        use MPI
        integer(kind=4) :: ia, ib, ierror
        real(kind=4) :: fa, fb
        call MPI_Reduce(ia, ib, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        call MPI_Reduce(fa, fb, 1, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        end"
    )
    check_fortran_source_compiles(${mpi_mismatched_args_prog} FORTRAN_MISMATCHED_MPI SRC_EXT F90)
    if(NOT FORTRAN_MISMATCHED_MPI)
        list(APPEND CMAKE_REQUIRED_FLAGS "-fallow-argument-mismatch")
        unset(FORTRAN_MISMATCHED_MPI CACHE)
        check_fortran_source_compiles(${mpi_mismatched_args_prog} FORTRAN_MISMATCHED_MPI SRC_EXT F90)

        if (FORTRAN_MISMATCHED_MPI)
            add_flags_if_avail(FLAGS "-fallow-argument-mismatch")
        endif()
    endif()

    # restore state
    if(DEFINED TMP_CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_FLAGS ${TMP_CMAKE_REQUIRED_FLAGS})
        unset(TMP_CMAKE_REQUIRED_FLAGS)
    else()
        unset(CMAKE_REQUIRED_FLAGS)
    endif()
    if(DEFINED TMP_CMAKE_REQUIRED_LIBRARIES)
        set(CMAKE_REQUIRED_LIBRARIES ${TMP_CMAKE_REQUIRED_LIBRARIES})
        unset(TMP_CMAKE_REQUIRED_LIBRARIES)
    else()
        unset(CMAKE_REQUIRED_LIBRARIES)
    endif()

    if(NOT FORTRAN_MISMATCHED_MPI)
        message(FATAL_ERROR "cannot compile Fortran calling MPI multiple times with different argument types")
    endif()

endfunction(allow_mismatched_mpi_args)
