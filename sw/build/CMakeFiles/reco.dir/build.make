# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build

# Include any dependencies generated for this target.
include CMakeFiles/reco.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/reco.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/reco.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/reco.dir/flags.make

CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o: CMakeFiles/reco.dir/flags.make
CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o: ../plugin/dRICH_reco.cxx
CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o: CMakeFiles/reco.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o -MF CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o.d -o CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o -c /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/dRICH_reco.cxx

CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/dRICH_reco.cxx > CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.i

CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/dRICH_reco.cxx -o CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.s

# Object files for target reco
reco_OBJECTS = \
"CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o"

# External object files for target reco
reco_EXTERNAL_OBJECTS =

reco: CMakeFiles/reco.dir/plugin/dRICH_reco.cxx.o
reco: CMakeFiles/reco.dir/build.make
reco: libdRICH_lib.a
reco: /usr/lib64/root/libCore.so
reco: /usr/lib64/root/libImt.so
reco: /usr/lib64/root/libRIO.so
reco: /usr/lib64/root/libNet.so
reco: /usr/lib64/root/libHist.so
reco: /usr/lib64/root/libGraf.so
reco: /usr/lib64/root/libGraf3d.so
reco: /usr/lib64/root/libGpad.so
reco: /usr/lib64/root/libROOTDataFrame.so
reco: /usr/lib64/root/libTree.so
reco: /usr/lib64/root/libTreePlayer.so
reco: /usr/lib64/root/libRint.so
reco: /usr/lib64/root/libPostscript.so
reco: /usr/lib64/root/libMatrix.so
reco: /usr/lib64/root/libPhysics.so
reco: /usr/lib64/root/libMathCore.so
reco: /usr/lib64/root/libThread.so
reco: /usr/lib64/root/libMultiProc.so
reco: /usr/lib64/root/libROOTVecOps.so
reco: CMakeFiles/reco.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable reco"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/reco.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/reco.dir/build: reco
.PHONY : CMakeFiles/reco.dir/build

CMakeFiles/reco.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/reco.dir/cmake_clean.cmake
.PHONY : CMakeFiles/reco.dir/clean

CMakeFiles/reco.dir/depend:
	cd /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles/reco.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/reco.dir/depend

