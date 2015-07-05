

set(command "/usr/local/Cellar/cmake/3.2.3/bin/cmake;-Dmake=${make};-Dconfig=${config};-P;/Users/tom/Code/KronMat/build/third_party/src/googletest-stamp/googletest-download-impl.cmake")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_FILE "/Users/tom/Code/KronMat/build/third_party/src/googletest-stamp/googletest-download-out.log"
  ERROR_FILE "/Users/tom/Code/KronMat/build/third_party/src/googletest-stamp/googletest-download-err.log"
  )
if(result)
  set(msg "Command failed: ${result}\n")
  foreach(arg IN LISTS command)
    set(msg "${msg} '${arg}'")
  endforeach()
  set(msg "${msg}\nSee also\n  /Users/tom/Code/KronMat/build/third_party/src/googletest-stamp/googletest-download-*.log")
  message(FATAL_ERROR "${msg}")
else()
  set(msg "googletest download command succeeded.  See also /Users/tom/Code/KronMat/build/third_party/src/googletest-stamp/googletest-download-*.log")
  message(STATUS "${msg}")
endif()
