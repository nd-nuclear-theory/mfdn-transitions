cmake_minimum_required(VERSION 3.20)

# Determine if mfdn-transitions is built as a subproject (using
# add_subdirectory) or if it is the master project.
if(NOT DEFINED MFDN_TRANSITIONS_MASTER_PROJECT)
  set(MFDN_TRANSITIONS_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(MFDN_TRANSITIONS_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

project(mfdn-transitions Fortran)

# make sure that the default is a Release
if(MFDN_TRANSITIONS_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel"
        FORCE
  )
endif()

# set location for Fortran module (.mod) files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

# git revision
find_package(Git REQUIRED QUIET)
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --long --tags --always --dirty
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE PACKAGE_GIT_VERSION
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
)
# add_compile_definitions(VCS_REVISION="${PACKAGE_GIT_VERSION}")
set_property(DIRECTORY PROPERTY CMAKE_CONFIGURE_DEPENDS .git/index)

# ##############################################################################
# configurable options
# ##############################################################################

set(TRANSITIONS_MAX_NUM_TBO
    8
    CACHE STRING "maximum number of two-body operators"
)
set(TRANSITIONS_MAX_NUM_KETS
    8
    CACHE STRING "maximum number of kets"
)
set(TRANSITIONS_MAX_PATHLENGTH
    1024
    CACHE STRING "maximum path length"
)
option(USE_AVX512 "Enable AVX-512 architecture flags" OFF)
option(USE_SANTITIZER "Enable leak detection" OFF)
option(USE_CHECK "Enable bounds checking" OFF)
option(USE_IPO "Enable interprocedural optimization/link-time optimization" ON)
option(GEN_OPT_REPORT "Generate optimization report" ON)

# ##############################################################################
# compiler flags/options
# ##############################################################################

include(cmake/compiler_flags.cmake)

# define preprocessor macros
add_compile_definitions(
  MaxNumTBMEops=${TRANSITIONS_MAX_NUM_TBO}
  MaxNumKets=${TRANSITIONS_MAX_NUM_KETS}
  MAX_PATHLENGTH=${TRANSITIONS_MAX_PATHLENGTH}
)

# enable preprocessor
set(CMAKE_Fortran_PREPROCESS ON)

# general optimization flags
add_flags_if_avail(FLAGS "-O3" CONFIGURATIONS Release RelWithDebInfo)
add_flags_if_avail(FLAGS "-O0" CONFIGURATIONS Debug)

# general debugging flags
add_flags_if_avail(FLAGS "-g" CONFIGURATIONS Debug RelWithDebInfo)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # debugging flags supported by ifort
  add_flags_if_avail(
    FLAGS "-warn" "-traceback"
    CONFIGURATIONS Debug RelWithDebInfo
  )
else()
  # optimization flags supported by GCC
  add_flags_if_avail(
    FLAGS "-ffast-math" "-funroll-loops"
    CONFIGURATIONS Release RelWithDebInfo
  )
  # debugging flags supported by GCC
  add_flags_if_avail(
    FLAGS "-W" "-Wall" "-Wextra" "-fbacktrace"
    CONFIGURATIONS Debug RelWithDebInfo
  )
endif()

# enable interprocedural optimization/link-time optimization
if(USE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(
    RESULT MFDN_TRANSITIONS_IPO_SUPPORTED
    OUTPUT MFDN_TRANSITIONS_IPO_SUPPORTED_MSG
  )
  if(MFDN_TRANSITIONS_IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
  else()
    message(
      WARNING "IPO is not supported: ${MFDN_TRANSITIONS_IPO_SUPPORTED_MSG}"
    )
  endif()
endif()

# enable runtime leak detection
if(USE_SANITIZER)
  add_flags_if_avail(FLAGS "-fsanitize=leak")

  # emit warnings that this feature is not available
  if(NOT Fortran_fsanitize_leak)
    message(WARNING "Sanitizer is not available for selected Fortran compiler")
  endif()
endif()

# check for AVX-512 flags
if(USE_AVX512)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_flags_if_avail(FLAGS "-xMIC-AVX512")
  else()
    # these flags are supported by newer GCC versions
    add_flags_if_avail(FLAGS "-mavx512f" "-mavx512pf" "-mavx512er" "-mavx512cd")
  endif()
endif()

# check for bounds-checking flags
if(USE_CHECK)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_flags_if_avail(FLAGS "-check all" CONFIGURATIONS Debug RelWithDebInfo)
  else()
    add_flags_if_avail(FLAGS "-fcheck=all" CONFIGURATIONS Debug RelWithDebInfo)
  endif()
endif()

# enable optimization report
if(GEN_OPT_REPORT)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_flags_if_avail(
      FLAGS "-list" "-qopt-report=1" "-qopt-report-phase=loop,openmp,vec"
    )
  else()
    add_flags_if_avail(FLAGS "-fopt-info-loop-omp-vec-all")
  endif()
endif()

# ##############################################################################
# find external projects/dependencies
# ##############################################################################

find_package(GSL REQUIRED)

include(cmake/allow_mismatched_mpi_args.cmake)

if(NOT TARGET MPI::MPI_Fortran)
  find_package(MPI REQUIRED)
endif()
allow_mismatched_mpi_args()

if(NOT TARGET OpenMP::OpenMP_Fortran)
  find_package(OpenMP REQUIRED)
endif()

# C needs to be enabled to find MKL BLAS
enable_language(C)
if(NOT TARGET BLAS::BLAS)
  find_package(BLAS REQUIRED)
endif()

# ##############################################################################
# define mfdn-transitions sources
# ##############################################################################

# define list of sources
set(mfdn-transitions_SOURCES
    src/module_Transitions_input.f90
    src/module_nodeinfo.f90
    src/module_SPbasis.f90
    src/module_Sparsity.f90
    src/module_Wig3J.f90
    src/module_MjStates.f90
    src/module_Wavefunctions.f90
    # src/module_TBME_Tz0.f90  # added only to mfdn-transitions target below
    # src/module_TBME_Tz12.f90 # added only to mfdn-transitions-DeltaTz target below
    src/module_TBops.f90
    src/module_OBops.f90
    src/subrts_general.f90
    src/subrts_bitrep.f90
    src/subrts_read_MBbasis.f90
    src/subrts_MBnonzero.f90
    src/subrts_OBDME.f90
    src/subrts_EvalOBDME.f90
    src/subrts_inputTBME.f90
    src/subrts_EvalTBTrans.f90
    src/subrts_ApplyTBTrans.f90
)

# define mfdn-transitions object library and add source files
add_library(mfdn-transitions OBJECT ${mfdn-transitions_SOURCES} src/module_TBME_Tz0.f90)
add_library(mfdn-transitions::mfdn-transitions ALIAS mfdn-transitions)

add_library(mfdn-transitions-DeltaTz OBJECT ${mfdn-transitions_SOURCES} src/module_TBME_Tz12.f90)
target_compile_definitions(mfdn-transitions-DeltaTz PUBLIC DeltaTz)
add_library(mfdn-transitions::mfdn-transitions-DeltaTz ALIAS mfdn-transitions-DeltaTz)

# ##############################################################################
# link dependencies and flags
# ##############################################################################

target_link_libraries(mfdn-transitions PUBLIC compile-options)
target_link_libraries(
  mfdn-transitions PUBLIC MPI::MPI_Fortran OpenMP::OpenMP_Fortran BLAS::BLAS GSL::gsl
)

target_link_libraries(mfdn-transitions-DeltaTz PUBLIC compile-options GSL::gsl)
target_link_libraries(
  mfdn-transitions-DeltaTz PUBLIC MPI::MPI_Fortran OpenMP::OpenMP_Fortran BLAS::BLAS GSL::gsl
)

# ##############################################################################
# define main executable
# ##############################################################################

# xtransitions
add_executable(xtransitions src/MFDn_Transitions.f90)
target_link_libraries(xtransitions PRIVATE compile-options)
target_link_libraries(
  xtransitions PRIVATE mfdn-transitions::mfdn-transitions MPI::MPI_Fortran
                       OpenMP::OpenMP_Fortran
)
target_compile_definitions(xtransitions PRIVATE VCS_REVISION="${PACKAGE_GIT_VERSION}")

# xtransitions-DeltaTz
add_executable(xtransitions-DeltaTz src/MFDn_Transitions.f90)
target_compile_definitions(xtransitions-DeltaTz PRIVATE DeltaTz)
target_link_libraries(xtransitions-DeltaTz PRIVATE compile-options)
target_link_libraries(
  xtransitions-DeltaTz PRIVATE mfdn-transitions::mfdn-transitions-DeltaTz
                       MPI::MPI_Fortran
                       OpenMP::OpenMP_Fortran
)
target_compile_definitions(xtransitions-DeltaTz PRIVATE VCS_REVISION="${PACKAGE_GIT_VERSION}")

# xapply
add_executable(xapply src/MFDn_apply_operator.f90)
target_link_libraries(xapply PRIVATE compile-options)
target_link_libraries(
  xapply PRIVATE mfdn-transitions::mfdn-transitions MPI::MPI_Fortran
                       OpenMP::OpenMP_Fortran
)
target_compile_definitions(xapply PRIVATE VCS_REVISION="${PACKAGE_GIT_VERSION}")

# xapply-DeltaTz
add_executable(xapply-DeltaTz src/MFDn_apply_operator.f90)
target_compile_definitions(xapply-DeltaTz PRIVATE DeltaTz)
target_link_libraries(xapply-DeltaTz PRIVATE compile-options)
target_link_libraries(
  xapply-DeltaTz PRIVATE mfdn-transitions::mfdn-transitions-DeltaTz
                       MPI::MPI_Fortran
                       OpenMP::OpenMP_Fortran
)
target_compile_definitions(xapply-DeltaTz PRIVATE VCS_REVISION="${PACKAGE_GIT_VERSION}")

# ##############################################################################
# define installation rules
#
# installation of targets is currently disabled, since the only exported target
# is an executable; if mfdn-transitions ever exports a library, this should be
# uncommented and used
# ##############################################################################

install(TARGETS xtransitions xtransitions-DeltaTz)
install(TARGETS xapply xapply-DeltaTz)

## install(TARGETS xtransitions EXPORT mfdn-transitionsTargets)
##
## install(
##   EXPORT mfdn-transitionsTargets
##   NAMESPACE mfdn-transitions::
##   FILE mfdn-transitionsTargets.cmake
##   DESTINATION lib/cmake/mfdn-transitions
## )
##
## include(CMakePackageConfigHelpers)
## # generate the config file that is includes the exports
## configure_package_config_file(
##   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
##   "${CMAKE_CURRENT_BINARY_DIR}/mfdn-transitionsConfig.cmake"
##   INSTALL_DESTINATION "lib/cmake/mfdn-transitions"
##   NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
## )
## # generate the version file for the config file
## write_basic_package_version_file(
##   "${CMAKE_CURRENT_BINARY_DIR}/mfdn-transitionsConfigVersion.cmake"
##   VERSION "${mfdn-transitions_VERSION_MAJOR}.${mfdn-transitions_VERSION_MINOR}"
##   COMPATIBILITY AnyNewerVersion
## )
##
## # install the configuration file
## install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mfdn-transitionsConfig.cmake
##               ${CMAKE_CURRENT_BINARY_DIR}/mfdn-transitionsConfigVersion.cmake
##         DESTINATION lib/cmake/mfdn-transitions
## )
##
## export(
##   EXPORT mfdn-transitionsTargets
##   NAMESPACE mfdn-transitions::
##   FILE "${CMAKE_CURRENT_BINARY_DIR}/mfdn-transitionsTargets.cmake"
## )

# ##############################################################################
# define tests
#
# tests are currently disabled, since no test programs are currently included
# in mfdn-transitions; if tests are implemented, this should be uncommented and
# used
# ##############################################################################

## set(mfdn-transitions_TEST_NAMES)
##
## add_custom_target(mfdn-transitions_tests)
## foreach(test_name IN LISTS mfdn-transitions_TEST_NAMES)
##   add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.f90)
##   target_link_libraries(${test_name} mfdn-transitions::mfdn-transitions)
##   add_dependencies(mfdn-transitions_tests ${test_name})
## endforeach()
##
## if(MFDN_TRANSITIONS_MASTER_PROJECT)
##   add_custom_target(tests)
##   add_dependencies(tests mfdn-transitions_tests)
## endif()
