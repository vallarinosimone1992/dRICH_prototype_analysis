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
include CMakeFiles/event.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/event.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/event.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/event.dir/flags.make

CMakeFiles/event.dir/plugin/event.cxx.o: CMakeFiles/event.dir/flags.make
CMakeFiles/event.dir/plugin/event.cxx.o: ../plugin/event.cxx
CMakeFiles/event.dir/plugin/event.cxx.o: CMakeFiles/event.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/event.dir/plugin/event.cxx.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/event.dir/plugin/event.cxx.o -MF CMakeFiles/event.dir/plugin/event.cxx.o.d -o CMakeFiles/event.dir/plugin/event.cxx.o -c /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/event.cxx

CMakeFiles/event.dir/plugin/event.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/event.dir/plugin/event.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/event.cxx > CMakeFiles/event.dir/plugin/event.cxx.i

CMakeFiles/event.dir/plugin/event.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/event.dir/plugin/event.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/plugin/event.cxx -o CMakeFiles/event.dir/plugin/event.cxx.s

# Object files for target event
event_OBJECTS = \
"CMakeFiles/event.dir/plugin/event.cxx.o"

# External object files for target event
event_EXTERNAL_OBJECTS =

event: CMakeFiles/event.dir/plugin/event.cxx.o
event: CMakeFiles/event.dir/build.make
event: libdRICH_lib.a
event: /usr/lib64/root/libCore.so
event: /usr/lib64/root/libImt.so
event: /usr/lib64/root/libRIO.so
event: /usr/lib64/root/libNet.so
event: /usr/lib64/root/libHist.so
event: /usr/lib64/root/libGraf.so
event: /usr/lib64/root/libGraf3d.so
event: /usr/lib64/root/libGpad.so
event: /usr/lib64/root/libROOTDataFrame.so
event: /usr/lib64/root/libTree.so
event: /usr/lib64/root/libTreePlayer.so
event: /usr/lib64/root/libRint.so
event: /usr/lib64/root/libPostscript.so
event: /usr/lib64/root/libMatrix.so
event: /usr/lib64/root/libPhysics.so
event: /usr/lib64/root/libMathCore.so
event: /usr/lib64/root/libThread.so
event: /usr/lib64/root/libMultiProc.so
event: /usr/lib64/root/libROOTVecOps.so
event: CMakeFiles/event.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable event"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/event.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/event.dir/build: event
.PHONY : CMakeFiles/event.dir/build

CMakeFiles/event.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/event.dir/cmake_clean.cmake
.PHONY : CMakeFiles/event.dir/clean

CMakeFiles/event.dir/depend:
	cd /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build /home/simone/Work/EIC/dRICH/prototype/testBeam/dRICH_prototype_analysis/sw/build/CMakeFiles/event.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/event.dir/depend

