set(dir ${CMAKE_CURRENT_BINARY_DIR}/test/${t})
execute_process(
  WORKING_DIRECTORY ${dir}
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/../test_problems/${in} ${in}
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../test_problems/truth/${t} truth
  COMMAND mpirun -n ${n} ${CMAKE_CURRENT_BINARY_DIR}/waveqlab3d ${in})
set(results 0)
execute_process(
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../python/read_binary.py ${dir} ${prefix}
  RESULT_VARIABLE test_fail)
math(EXPR results "${results} + ${test_fail}")
if (results)
  message( SEND_ERROR "waveqlab3d output for ${in} does not match true solution." )
endif (results)