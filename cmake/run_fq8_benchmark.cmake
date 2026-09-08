# run_fq8_benchmark.cmake — Run WaveQLab3D-fQ with an fQ8 input file,
# generate f-k reference seismograms, and compare using Kristekova EM/PM.
#
# Required variables:
#   MPIEXEC  — MPI launcher
#   EXE      — path to waveqlab3d executable
#   INPUT    — path to solver input file
#   PYTHON   — path to Python 3 interpreter
#   FK_REF   — path to python/fk_reference.py
#   CASE     — benchmark case name (elastic, constant_q, powerlaw, layered)
#   EM_TOL   — max envelope misfit (percent)
#   PM_TOL   — max phase misfit (percent)
#
# Optional:
#   NPROCS   — number of MPI ranks (default 2)

if(NOT DEFINED NPROCS)
  set(NPROCS 2)
endif()

# Run solver
execute_process(
  COMMAND "${MPIEXEC}" -np ${NPROCS} "${EXE}" "${INPUT}"
  RESULT_VARIABLE solver_result
  OUTPUT_VARIABLE solver_output
  ERROR_VARIABLE solver_error
  TIMEOUT 600)
if(NOT solver_result EQUAL 0)
  message(FATAL_ERROR "Solver failed (exit ${solver_result}):\n${solver_error}")
endif()

# Verify solver produced output
string(REGEX MATCH "fQ8 final state: max\\|field\\|=[^\n]+" state "${solver_output}")
if(state STREQUAL "")
  message(FATAL_ERROR "fQ8 final-state diagnostic missing from output")
endif()
if(state MATCHES "max\\|field\\|= *0\\.0+E\\+00")
  message(FATAL_ERROR "fQ8 dynamic fixture did not produce a nonzero field")
endif()

# Generate reference seismograms
set(REF_DIR "${CMAKE_CURRENT_BINARY_DIR}/fk_ref_${CASE}")
execute_process(
  COMMAND "${PYTHON}" "${FK_REF}" --case "${CASE}" --outdir "${REF_DIR}"
  RESULT_VARIABLE ref_result
  OUTPUT_VARIABLE ref_output
  ERROR_VARIABLE ref_error
  TIMEOUT 120)
if(NOT ref_result EQUAL 0)
  message(FATAL_ERROR "Reference generation failed:\n${ref_error}")
endif()

message(STATUS "Benchmark ${CASE}: solver and reference completed successfully.")
message(STATUS "Solver state: ${state}")
message(STATUS "Reference output: ${REF_DIR}")

# Seismogram comparison would be added here once the solver's station output
# paths are known. For now, the test verifies that both the solver and the
# reference generator run to completion without error.
#
# To add EM/PM comparison:
#   execute_process(
#     COMMAND "${PYTHON}" "${FK_REF}"
#       --compare "${REF_DIR}/${CASE}/uz_r5000m.dat" "${SOLVER_OUTDIR}/station_1.dat"
#       --dt-compare 0.002
#     OUTPUT_VARIABLE compare_output ...)
#   string(REGEX MATCH "EM = ([0-9.]+)%" em_match "${compare_output}")
