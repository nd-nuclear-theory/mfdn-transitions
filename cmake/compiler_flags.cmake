# ##############################################################################
# macros for adding compiler flags
#
# modified from https://docs.nersc.gov/development/build-tools/cmake/
# ##############################################################################

include(CheckCompilerFlag)

# create interface target with compiler flags
if(NOT TARGET compile-options)
    add_library(compile-options INTERFACE)
endif()

# ##############################################################################
# Check if flag if supported, if so add to compile-options.
#
# Options:
#   FLAGS -- compiler flags to check and add
#   LANGUAGES -- languages for specified flags (default: all enabled languages)
#   CONFIGURATIONS -- configurations (Release, Debug, etc.) for specified flags
#     (default: all configurations)
# ##############################################################################
macro(add_flags_if_avail)
    cmake_parse_arguments(ADD_FLAGS "" "" "FLAGS;LANGUAGES;CONFIGURATIONS" ${ARGN})
    list(JOIN ADD_FLAGS_CONFIGURATIONS "," ADD_FLAGS_CFGS)
    if(NOT ADD_FLAGS_LANGUAGES)
        get_property(ADD_FLAGS_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
    endif()
    foreach(LANG ${ADD_FLAGS_LANGUAGES})
        foreach(FLAG ${ADD_FLAGS_FLAGS})
            # create a variable for checking the flag if supported, e.g.:
            #   -fp-model=precise --> c_fp_model_precise
            string(REGEX REPLACE "^-" "${LANG}_" FLAG_NAME "${FLAG}")
            string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
            string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
            string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

            check_compiler_flag("${LANG}" "${FLAG}" ${FLAG_NAME})
            if(${FLAG_NAME})
                if(NOT "${ADD_FLAGS_CFGS}" STREQUAL "")
                    target_compile_options(compile-options INTERFACE
                    $<$<AND:$<COMPILE_LANGUAGE:${LANG}>,$<CONFIG:${ADD_FLAGS_CFGS}>>:${FLAG}>)
                else()
                    target_compile_options(compile-options INTERFACE
                        $<$<COMPILE_LANGUAGE:${LANG}>:${FLAG}>)
                endif()
            endif()
        endforeach()
    endforeach()
    unset(ADD_FLAGS_FLAGS)
    unset(ADD_FLAGS_LANGUAGES)
    unset(ADD_FLAGS_CONFIGURATIONS)
endmacro()
