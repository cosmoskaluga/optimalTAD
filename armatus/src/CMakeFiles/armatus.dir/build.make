# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/anna/Downloads/tads/ArmatusParallel/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/anna/Downloads/tads/ArmatusParallel/src

# Include any dependencies generated for this target.
include CMakeFiles/armatus.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/armatus.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/armatus.dir/flags.make

CMakeFiles/armatus.dir/Armatus.cpp.o: CMakeFiles/armatus.dir/flags.make
CMakeFiles/armatus.dir/Armatus.cpp.o: Armatus.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/armatus.dir/Armatus.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/armatus.dir/Armatus.cpp.o -c /home/anna/Downloads/tads/ArmatusParallel/src/Armatus.cpp

CMakeFiles/armatus.dir/Armatus.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armatus.dir/Armatus.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anna/Downloads/tads/ArmatusParallel/src/Armatus.cpp > CMakeFiles/armatus.dir/Armatus.cpp.i

CMakeFiles/armatus.dir/Armatus.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armatus.dir/Armatus.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anna/Downloads/tads/ArmatusParallel/src/Armatus.cpp -o CMakeFiles/armatus.dir/Armatus.cpp.s

CMakeFiles/armatus.dir/ArmatusUtil.cpp.o: CMakeFiles/armatus.dir/flags.make
CMakeFiles/armatus.dir/ArmatusUtil.cpp.o: ArmatusUtil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/armatus.dir/ArmatusUtil.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/armatus.dir/ArmatusUtil.cpp.o -c /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusUtil.cpp

CMakeFiles/armatus.dir/ArmatusUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armatus.dir/ArmatusUtil.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusUtil.cpp > CMakeFiles/armatus.dir/ArmatusUtil.cpp.i

CMakeFiles/armatus.dir/ArmatusUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armatus.dir/ArmatusUtil.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusUtil.cpp -o CMakeFiles/armatus.dir/ArmatusUtil.cpp.s

CMakeFiles/armatus.dir/ArmatusParams.cpp.o: CMakeFiles/armatus.dir/flags.make
CMakeFiles/armatus.dir/ArmatusParams.cpp.o: ArmatusParams.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/armatus.dir/ArmatusParams.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/armatus.dir/ArmatusParams.cpp.o -c /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusParams.cpp

CMakeFiles/armatus.dir/ArmatusParams.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armatus.dir/ArmatusParams.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusParams.cpp > CMakeFiles/armatus.dir/ArmatusParams.cpp.i

CMakeFiles/armatus.dir/ArmatusParams.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armatus.dir/ArmatusParams.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusParams.cpp -o CMakeFiles/armatus.dir/ArmatusParams.cpp.s

CMakeFiles/armatus.dir/ArmatusDAG.cpp.o: CMakeFiles/armatus.dir/flags.make
CMakeFiles/armatus.dir/ArmatusDAG.cpp.o: ArmatusDAG.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/armatus.dir/ArmatusDAG.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/armatus.dir/ArmatusDAG.cpp.o -c /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusDAG.cpp

CMakeFiles/armatus.dir/ArmatusDAG.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armatus.dir/ArmatusDAG.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusDAG.cpp > CMakeFiles/armatus.dir/ArmatusDAG.cpp.i

CMakeFiles/armatus.dir/ArmatusDAG.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armatus.dir/ArmatusDAG.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anna/Downloads/tads/ArmatusParallel/src/ArmatusDAG.cpp -o CMakeFiles/armatus.dir/ArmatusDAG.cpp.s

CMakeFiles/armatus.dir/IntervalScheduling.cpp.o: CMakeFiles/armatus.dir/flags.make
CMakeFiles/armatus.dir/IntervalScheduling.cpp.o: IntervalScheduling.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/armatus.dir/IntervalScheduling.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/armatus.dir/IntervalScheduling.cpp.o -c /home/anna/Downloads/tads/ArmatusParallel/src/IntervalScheduling.cpp

CMakeFiles/armatus.dir/IntervalScheduling.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armatus.dir/IntervalScheduling.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anna/Downloads/tads/ArmatusParallel/src/IntervalScheduling.cpp > CMakeFiles/armatus.dir/IntervalScheduling.cpp.i

CMakeFiles/armatus.dir/IntervalScheduling.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armatus.dir/IntervalScheduling.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anna/Downloads/tads/ArmatusParallel/src/IntervalScheduling.cpp -o CMakeFiles/armatus.dir/IntervalScheduling.cpp.s

# Object files for target armatus
armatus_OBJECTS = \
"CMakeFiles/armatus.dir/Armatus.cpp.o" \
"CMakeFiles/armatus.dir/ArmatusUtil.cpp.o" \
"CMakeFiles/armatus.dir/ArmatusParams.cpp.o" \
"CMakeFiles/armatus.dir/ArmatusDAG.cpp.o" \
"CMakeFiles/armatus.dir/IntervalScheduling.cpp.o"

# External object files for target armatus
armatus_EXTERNAL_OBJECTS =

armatus: CMakeFiles/armatus.dir/Armatus.cpp.o
armatus: CMakeFiles/armatus.dir/ArmatusUtil.cpp.o
armatus: CMakeFiles/armatus.dir/ArmatusParams.cpp.o
armatus: CMakeFiles/armatus.dir/ArmatusDAG.cpp.o
armatus: CMakeFiles/armatus.dir/IntervalScheduling.cpp.o
armatus: CMakeFiles/armatus.dir/build.make
armatus: CMakeFiles/armatus.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable armatus"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/armatus.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/armatus.dir/build: armatus

.PHONY : CMakeFiles/armatus.dir/build

CMakeFiles/armatus.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/armatus.dir/cmake_clean.cmake
.PHONY : CMakeFiles/armatus.dir/clean

CMakeFiles/armatus.dir/depend:
	cd /home/anna/Downloads/tads/ArmatusParallel/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/anna/Downloads/tads/ArmatusParallel/src /home/anna/Downloads/tads/ArmatusParallel/src /home/anna/Downloads/tads/ArmatusParallel/src /home/anna/Downloads/tads/ArmatusParallel/src /home/anna/Downloads/tads/ArmatusParallel/src/CMakeFiles/armatus.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/armatus.dir/depend

