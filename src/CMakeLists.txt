cmake_minimum_required(VERSION 3.5)
project (WaveQLab3D Fortran)
set(CMAKE_VERBOSE_MAKEFILE ON)

# This branch requires MPI for compiling
FIND_PACKAGE(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARY_PATH})


# This branch requires HDF5
#FIND_PACKAGE(HDF5 COMPONENTS Fortran)
#add_definitions(${HDF5_DEFINITIONS})
#include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# The fortran compiler flags depends on the compiler being used
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-fdefault-real-8 -fdefault-double-8 -O5 -Wuninitialized")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fdefault-real-8 -fdefault-double-8 -Og -g -Wall -Wextra -Wconversion -Wuninitialized -fbounds-check -fbacktrace -finit-real=nan -ffpe-trap=invalid,zero,overflow -fimplicit-none -fall-intrinsics -std=f2008ts")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-m64 -xHost -r8 -i4 -std08 -g -O2 -heap-arrays")
    # Keep runtime bounds/interface/stack checks usable.  -ftrapuv together
    # with -fpe0 aborts before main while legacy derived types are constructed.
    set (CMAKE_Fortran_FLAGS_DEBUG   "-r8 -i4 -std08 -g -heap-arrays -check all,noarg_temp_created -warn all -traceback -fp-stack-check -fltconsistency")
else (Fortran_COMPILER_NAME MATCHES "ftn.*")
    # Cray fortran compiler
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -Oipa0 -v -V")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# This needs to be set for the source code to become visible in the debugger.
# It will need to be expanded when we start using Linux and other OSes.
if(APPLE)
    set(DEBUGGING_SYMBOLS YES)
    set(DEBUG_INFORMATION_FORMAT dwarf-with-dsym)
endif(APPLE)

# Setup the source files and the executable name for the project.
set(WQL3D_SRC
  fault_output.f90
  fault_receiver_output.f90
  iface.f90
  block.f90
  boundary.f90
  common.f90
  diagnostics.f90
  decomposition_safety.f90
  datatypes.f90
  domain.f90
  metrics.f90
  elastic.f90
  plastic.f90
  fields.f90
  grid.f90
  main.f90
  withers_tables.f90
  anelastic_q4_model.f90
  anelastic_q8_model.f90
  anelastic_cq8_b2_model.f90
  anelastic_cq_model.f90
  anelastic_fq8_model.f90
  simulation_config.f90
  input_preflight.f90
  material.f90
  anelastic_cq_material.f90
  plastic_material.f90
  pml.f90
  mpi3dbasic.f90
  mpi3dcomm.f90
  mpi3dio.f90
  mpi3d_interface.f90
  time_step.f90
  JU_xJU_yJU_z6.f90
  RHS_Interior.f90
  mms.f90
  BoundaryConditions.f90
  CouplingForcing.f90
  Interface_Condition.f90
  initial_stress_condition.f90
  inter_material.f90
  unit_normals.f90
  seismogram.f90
  slice_output.f90
  plane_output.f90
  #hdf5_output.f90
  parallel_write.f90
  moment_tensor.f90
  )

set (PRE_WQL3D_SRC
    common.f90
    diagnostics.f90
    datatypes.f90
    grid.f90
    preprocessor.f90
    metrics.f90
    mpi3dbasic.f90
    mpi3dcomm.f90
    mpi3dio.f90
    mpi3d_interface.f90
  )

add_executable (waveqlab3d ${WQL3D_SRC})
add_executable (pre_wql3d ${PRE_WQL3D_SRC})
target_link_libraries (waveqlab3d ${MPI_Fortran_LIBRARIES})
target_link_libraries (pre_wql3d ${MPI_Fortran_LIBRARIES})
#target_link_libraries (waveqlab3d ${HDF5_LIBRARIES})

################################################
# Enable testing using example files
###############################################

enable_testing()

add_executable(fq8_effective_response_test
  common.f90
  withers_tables.f90
  anelastic_fq8_model.f90
  ../tests/fq8_effective_response_test.f90)
set_target_properties(fq8_effective_response_test PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fq8_test_modules)
target_include_directories(fq8_effective_response_test PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/fq8_test_modules)
add_test(NAME fq8_effective_response_unit
  COMMAND $<TARGET_FILE:fq8_effective_response_test>)

add_executable(cq_coefficients_test
  common.f90
  anelastic_cq_model.f90
  ../tests/cq_coefficients_test.f90)
set_target_properties(cq_coefficients_test PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cq_test_modules)
target_include_directories(cq_coefficients_test PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/cq_test_modules)
add_test(NAME cq_coefficients_unit COMMAND $<TARGET_FILE:cq_coefficients_test>)

# Intel MPI's mpirun wrapper supplies its runtime-library environment, while
# the sibling mpiexec binary on this cluster does not.
get_filename_component(MPIEXEC_DIRECTORY ${MPIEXEC_EXECUTABLE} DIRECTORY)
set(WQL3D_MPI_RUNNER ${MPIEXEC_DIRECTORY}/mpirun)


set (no_mpi_proc 4)

file(MAKE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/test/1/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/2/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/3/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/4/data)
  
add_test(NAME mpi_test1
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_cart_fric.in
    -D prefix=test_rup_cart_fric
    -D t=1
    -D n=${no_mpi_proc}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)

add_test(NAME q8_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_dynamic.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq8_b2_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ8_b2_dynamic.in
    -D STATE_LABEL=cQ8-b2
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq8_b2_rejects_one_block
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D NPROCS=1
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ8_b2_invalid_nb1.in
    -D EXPECTED_CODE=CFG-CQ8-B2-002
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

add_test(NAME cq_mechanisms_4_through_8
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_mechanism_regression.cmake)

add_test(NAME cq_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in
    -D STATE_LABEL=cQ
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq_rejects_invalid_mechanism_counts
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_invalid_regression.cmake)

add_test(NAME cq_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_pml_regression.cmake)

add_test(NAME cq_traditional_and_drp_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_fd_regression.cmake)

add_test(NAME q4_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_dynamic.in
    -D STATE_LABEL=Q4
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME q4_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_dynamic_pml.in
    -D STATE_LABEL=Q4
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME q4_asymmetric_three_ranks
  COMMAND ${WQL3D_MPI_RUNNER} ${MPIEXEC_NUMPROC_FLAG} 3
    $<TARGET_FILE:waveqlab3d>
    ${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_asymmetric.in)

add_test(NAME q4_deprecated_alias
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_alias.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q4_alias_regression.cmake)

foreach(q4_invalid_case legacy partial method nonpositive band)
  add_test(NAME q4_rejects_${q4_invalid_case}
    COMMAND ${CMAKE_COMMAND}
      -D MPIEXEC=${WQL3D_MPI_RUNNER}
      -D NPROCS=2
      -D EXE=$<TARGET_FILE:waveqlab3d>
      -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_invalid_${q4_invalid_case}.in
      -D EXPECTED_CODE=CFG-Q4-001
      -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)
endforeach()

add_test(NAME q8_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_dynamic_pml.in
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_coarse_grain_2_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_withers.in
    -D STATE_LABEL=fQ8
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_coarse_grain_0_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_full.in
    -D STATE_LABEL=fQ8
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_rejects_full_withers
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D NPROCS=2
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_invalid_full_withers.in
    -D EXPECTED_CODE=CFG-FQ8-001
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

add_test(NAME q8_asymmetric_three_ranks
  COMMAND ${WQL3D_MPI_RUNNER} ${MPIEXEC_NUMPROC_FLAG} 3
    $<TARGET_FILE:waveqlab3d>
    ${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_asymmetric.in)

add_test(NAME decomposition_rejects_excessive_split
  COMMAND ${CMAKE_COMMAND}
    -D MPIEXEC=${WQL3D_MPI_RUNNER}
    -D NPROCS=4
    -D EXE=$<TARGET_FILE:waveqlab3d>
    -D INPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_current.in
    -D EXPECTED_CODE=CFG-DECOMP-101
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

add_test(NAME mpi_test2
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_curv_fric.in
    -D prefix=test_rup_curv_fric
    -D t=2
    -D n=${no_mpi_proc}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)

add_test(NAME premesh_test1
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_cart_fric_premesh.in
    -D prefix=test_rup_cart_fric
    -D t=1
    -D n=${no_mpi_proc}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_test.cmake)

add_test(NAME premesh_test2
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_curv_fric_premesh.in
    -D prefix=test_rup_curv_fric
    -D t=2
    -D n=${no_mpi_proc}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_test.cmake)

add_test(NAME serial_test1
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_cart_fric.in
    -D prefix=test_rup_cart_fric
    -D t=1
    -D n=1
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)

add_test(NAME serial_test2
  COMMAND ${CMAKE_COMMAND} 
    -D in=test_rup_curv_fric.in
    -D prefix=test_rup_curv_fric
    -D t=2
    -D n=1
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)

##############################################################
#
# Add a target to generate API documentation with Doxygen
#
##############################################################

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../conf/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(docs
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)
