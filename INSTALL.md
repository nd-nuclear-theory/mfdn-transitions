Notes for installing mfdn-transitions with cmake

06/16/23 (mac): Extract from README.md and add note on cmake --install.

----------------------------------------------------------------

To build this code using CMake, first you must set up the cmake build directory:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cmake -B build/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then use CMake to carry out the build (i.e., to run the makefile):

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cmake --build build/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to use the postprocessor with the `mcscript-ncci` scripting, you now
need to install the executable to the location expected by the scripting.  You
will need to use the appropriate install prefix.  In a non-Cray environment:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cmake --install build/ --prefix=${MCSCRIPT_INSTALL_HOME}/mfdn-transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Or, in a Cray environment (such as NERSC), where we have a different installation
directory for each machine architecture:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cmake --install build/ --prefix=${MCSCRIPT_INSTALL_HOME}/${CRAY_CPU_TARGET}/mfdn-transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

