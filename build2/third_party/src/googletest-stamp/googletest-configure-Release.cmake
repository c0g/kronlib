

set(command "/usr/local/Cellar/cmake/3.3.2/bin/cmake;-GUnix Makefiles;/Users/tom/Code/KronMat/build2/third_party/src/googletest")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_FILE "/Users/tom/Code/KronMat/build2/third_party/src/googletest-stamp/googletest-configure-out.log"
  ERROR_FILE "/Users/tom/Code/KronMat/build2/third_party/src/googletest-stamp/googletest-configure-err.log"
  )
if(result)
  set(msg "Command failed: ${result}\n")
  foreach(arg IN LISTS command)
    set(msg "${msg} '${arg}'")
  endforeach()
  set(msg "${msg}\nSee also\n  /Users/tom/Code/KronMat/build2/third_party/src/googletest-stamp/googletest-configure-*.log")
  message(FATAL_ERROR "${msg}")
else()
  set(msg "googletest configure command succeeded.  See also /Users/tom/Code/KronMat/build2/third_party/src/googletest-stamp/googletest-configure-*.log")
  message(STATUS "${msg}")
endif()
