# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jmahler/Libraries/phat_1_4_0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jmahler/Libraries/phat_1_4_0/build

# Include any dependencies generated for this target.
include CMakeFiles/simple_example.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/simple_example.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simple_example.dir/flags.make

CMakeFiles/simple_example.dir/src/simple_example.cpp.o: CMakeFiles/simple_example.dir/flags.make
CMakeFiles/simple_example.dir/src/simple_example.cpp.o: ../src/simple_example.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jmahler/Libraries/phat_1_4_0/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simple_example.dir/src/simple_example.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simple_example.dir/src/simple_example.cpp.o -c /home/jmahler/Libraries/phat_1_4_0/src/simple_example.cpp

CMakeFiles/simple_example.dir/src/simple_example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simple_example.dir/src/simple_example.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jmahler/Libraries/phat_1_4_0/src/simple_example.cpp > CMakeFiles/simple_example.dir/src/simple_example.cpp.i

CMakeFiles/simple_example.dir/src/simple_example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simple_example.dir/src/simple_example.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jmahler/Libraries/phat_1_4_0/src/simple_example.cpp -o CMakeFiles/simple_example.dir/src/simple_example.cpp.s

CMakeFiles/simple_example.dir/src/simple_example.cpp.o.requires:
.PHONY : CMakeFiles/simple_example.dir/src/simple_example.cpp.o.requires

CMakeFiles/simple_example.dir/src/simple_example.cpp.o.provides: CMakeFiles/simple_example.dir/src/simple_example.cpp.o.requires
	$(MAKE) -f CMakeFiles/simple_example.dir/build.make CMakeFiles/simple_example.dir/src/simple_example.cpp.o.provides.build
.PHONY : CMakeFiles/simple_example.dir/src/simple_example.cpp.o.provides

CMakeFiles/simple_example.dir/src/simple_example.cpp.o.provides.build: CMakeFiles/simple_example.dir/src/simple_example.cpp.o

# Object files for target simple_example
simple_example_OBJECTS = \
"CMakeFiles/simple_example.dir/src/simple_example.cpp.o"

# External object files for target simple_example
simple_example_EXTERNAL_OBJECTS =

simple_example: CMakeFiles/simple_example.dir/src/simple_example.cpp.o
simple_example: CMakeFiles/simple_example.dir/build.make
simple_example: CMakeFiles/simple_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable simple_example"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simple_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simple_example.dir/build: simple_example
.PHONY : CMakeFiles/simple_example.dir/build

CMakeFiles/simple_example.dir/requires: CMakeFiles/simple_example.dir/src/simple_example.cpp.o.requires
.PHONY : CMakeFiles/simple_example.dir/requires

CMakeFiles/simple_example.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simple_example.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simple_example.dir/clean

CMakeFiles/simple_example.dir/depend:
	cd /home/jmahler/Libraries/phat_1_4_0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jmahler/Libraries/phat_1_4_0 /home/jmahler/Libraries/phat_1_4_0 /home/jmahler/Libraries/phat_1_4_0/build /home/jmahler/Libraries/phat_1_4_0/build /home/jmahler/Libraries/phat_1_4_0/build/CMakeFiles/simple_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simple_example.dir/depend

