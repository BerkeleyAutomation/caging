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
include CMakeFiles/phat.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/phat.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/phat.dir/flags.make

CMakeFiles/phat.dir/src/phat.cpp.o: CMakeFiles/phat.dir/flags.make
CMakeFiles/phat.dir/src/phat.cpp.o: ../src/phat.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jmahler/Libraries/phat_1_4_0/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/phat.dir/src/phat.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/phat.dir/src/phat.cpp.o -c /home/jmahler/Libraries/phat_1_4_0/src/phat.cpp

CMakeFiles/phat.dir/src/phat.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phat.dir/src/phat.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jmahler/Libraries/phat_1_4_0/src/phat.cpp > CMakeFiles/phat.dir/src/phat.cpp.i

CMakeFiles/phat.dir/src/phat.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phat.dir/src/phat.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jmahler/Libraries/phat_1_4_0/src/phat.cpp -o CMakeFiles/phat.dir/src/phat.cpp.s

CMakeFiles/phat.dir/src/phat.cpp.o.requires:
.PHONY : CMakeFiles/phat.dir/src/phat.cpp.o.requires

CMakeFiles/phat.dir/src/phat.cpp.o.provides: CMakeFiles/phat.dir/src/phat.cpp.o.requires
	$(MAKE) -f CMakeFiles/phat.dir/build.make CMakeFiles/phat.dir/src/phat.cpp.o.provides.build
.PHONY : CMakeFiles/phat.dir/src/phat.cpp.o.provides

CMakeFiles/phat.dir/src/phat.cpp.o.provides.build: CMakeFiles/phat.dir/src/phat.cpp.o

# Object files for target phat
phat_OBJECTS = \
"CMakeFiles/phat.dir/src/phat.cpp.o"

# External object files for target phat
phat_EXTERNAL_OBJECTS =

phat: CMakeFiles/phat.dir/src/phat.cpp.o
phat: CMakeFiles/phat.dir/build.make
phat: CMakeFiles/phat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable phat"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/phat.dir/build: phat
.PHONY : CMakeFiles/phat.dir/build

CMakeFiles/phat.dir/requires: CMakeFiles/phat.dir/src/phat.cpp.o.requires
.PHONY : CMakeFiles/phat.dir/requires

CMakeFiles/phat.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/phat.dir/cmake_clean.cmake
.PHONY : CMakeFiles/phat.dir/clean

CMakeFiles/phat.dir/depend:
	cd /home/jmahler/Libraries/phat_1_4_0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jmahler/Libraries/phat_1_4_0 /home/jmahler/Libraries/phat_1_4_0 /home/jmahler/Libraries/phat_1_4_0/build /home/jmahler/Libraries/phat_1_4_0/build /home/jmahler/Libraries/phat_1_4_0/build/CMakeFiles/phat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/phat.dir/depend

