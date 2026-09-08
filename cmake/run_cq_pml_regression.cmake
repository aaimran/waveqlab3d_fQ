file(READ "${INPUT}" input_text)
string(REPLACE "btp(1)%pml_lqrs=F,F,F" "btp(1)%pml_lqrs=T,F,F" pml_text "${input_text}")
string(REPLACE "btp(2)%pml_rqrs=F,F,F" "btp(2)%pml_rqrs=T,F,F" pml_text "${pml_text}")
string(REPLACE "%npml=0" "%npml=5" pml_text "${pml_text}")
set(pml_input "${CMAKE_CURRENT_BINARY_DIR}/test_anelastic_cQ_dynamic_pml.in")
file(WRITE "${pml_input}" "${pml_text}")
execute_process(COMMAND "${CMAKE_COMMAND}"
  -D MPIEXEC=${MPIEXEC} -D EXE=${EXE} -D INPUT=${pml_input} -D STATE_LABEL=cQ
  -P ${CMAKE_CURRENT_LIST_DIR}/run_q8_dynamic_regression.cmake
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "anelastic-cQ PML regression failed: ${error}\n${output}")
endif()
