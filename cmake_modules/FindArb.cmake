if(Arb_INCLUDE_DIR AND Arb_LIBRARIES)
    # Already in cache, be silent
    set(Arb_FIND_QUIETLY TRUE)
endif()

find_path(Arb_INCLUDE_DIR arb.h)
find_library(Arb_LIBRARIES NAMES arb)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Arb DEFAULT_MSG Arb_INCLUDE_DIR Arb_LIBRARIES)

mark_as_advanced(Arb_INCLUDE_DIR Arb_LIBRARIES)
