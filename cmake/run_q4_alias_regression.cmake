execute_process(
  COMMAND "${MPIEXEC}" -np 2 "${EXE}" "${INPUT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "deprecated Q4 alias run failed: ${error}")
endif()
string(REGEX MATCHALL "CFG-Q4-DEP-001" warnings "${output}")
list(LENGTH warnings warning_count)
if(NOT warning_count EQUAL 1)
  message(FATAL_ERROR "expected exactly one CFG-Q4-DEP-001 warning, found ${warning_count}")
endif()
if(NOT output MATCHES "anelastic-Q4 parameters:")
  message(FATAL_ERROR "deprecated alias was not normalized to anelastic-Q4")
endif()
