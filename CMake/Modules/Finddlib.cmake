find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  pkg_search_module(DLIB dlib)
else()

  find_path(DLIB_ROOT_DIR
      dlib/algs.h
      PATHS "$ENV{BuildEnvironment}/dlib/17.47/"
      NO_DEFAULT_PATH)

  set(DLIB_INCLUDE_DIRS "${DLIB_ROOT_DIR}")
  set(DLIB_LIBRARY_DIRS "${DLIB_ROOT_DIR}/lib")

  set(DLIB_LIBRARIES "${DLIB_LIBRARY_DIRS}/dlib.lib")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(DLIB DEFAULT_MSG DLIB_INCLUDE_DIRS DLIB_LIBRARY_DIRS)

endif()