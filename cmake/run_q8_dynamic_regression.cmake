if(NOT DEFINED STATE_LABEL)
  set(STATE_LABEL Q8)
endif()

execute_process(
  COMMAND "${MPIEXEC}" -np 1 "${EXE}" "${INPUT}"
  RESULT_VARIABLE result1
  OUTPUT_VARIABLE output1
  ERROR_VARIABLE error1)
if(NOT result1 EQUAL 0)
  message(FATAL_ERROR "one-rank ${STATE_LABEL} dynamic run failed: ${error1}")
endif()

execute_process(
  COMMAND "${MPIEXEC}" -np 2 "${EXE}" "${INPUT}"
  RESULT_VARIABLE result2
  OUTPUT_VARIABLE output2
  ERROR_VARIABLE error2)
if(NOT result2 EQUAL 0)
  message(FATAL_ERROR "two-rank ${STATE_LABEL} dynamic run failed: ${error2}")
endif()

string(REGEX MATCH "${STATE_LABEL} final state: max\\|field\\|=[^\n]+" state1 "${output1}")
string(REGEX MATCH "${STATE_LABEL} final state: max\\|field\\|=[^\n]+" state2 "${output2}")
if(state1 STREQUAL "" OR state2 STREQUAL "")
  message(FATAL_ERROR "${STATE_LABEL} final-state diagnostic is missing")
endif()
if(NOT state1 STREQUAL state2)
  message(FATAL_ERROR "decomposition mismatch:\n1 rank: ${state1}\n2 ranks: ${state2}")
endif()
if(state1 MATCHES "max\\|field\\|= *0\\.0+E\\+00")
  message(FATAL_ERROR "${STATE_LABEL} dynamic fixture did not produce a nonzero field")
endif()
if(state1 MATCHES "max\\|memory\\|= *0\\.0+E\\+00")
  message(FATAL_ERROR "${STATE_LABEL} dynamic fixture did not evolve memory variables")
endif()

if(STATE_LABEL STREQUAL "cQ8-b2")
  string(FIND "${output1}" "block 1: Qs0 =   4.0000E+01, Qp0 =   8.0000E+01" block1_q_position)
  string(FIND "${output1}" "block 2: Qs0 =   1.0000E+02, Qp0 =   1.8000E+02" block2_q_position)
  if(block1_q_position EQUAL -1 OR block2_q_position EQUAL -1)
    message(FATAL_ERROR "cQ8-b2 output does not contain the resolved per-block Q mapping")
  endif()
endif()
