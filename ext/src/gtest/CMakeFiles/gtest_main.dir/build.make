# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /ssd/ig_repertoire_constructor/ext/src/gtest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ssd/ig_repertoire_constructor/ext/src/gtest

# Include any dependencies generated for this target.
include CMakeFiles/gtest_main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/gtest_main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gtest_main.dir/flags.make

CMakeFiles/gtest_main.dir/gtest_main.cc.o: CMakeFiles/gtest_main.dir/flags.make
CMakeFiles/gtest_main.dir/gtest_main.cc.o: gtest_main.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /ssd/ig_repertoire_constructor/ext/src/gtest/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/gtest_main.dir/gtest_main.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/gtest_main.dir/gtest_main.cc.o -c /ssd/ig_repertoire_constructor/ext/src/gtest/gtest_main.cc

CMakeFiles/gtest_main.dir/gtest_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gtest_main.dir/gtest_main.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /ssd/ig_repertoire_constructor/ext/src/gtest/gtest_main.cc > CMakeFiles/gtest_main.dir/gtest_main.cc.i

CMakeFiles/gtest_main.dir/gtest_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gtest_main.dir/gtest_main.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /ssd/ig_repertoire_constructor/ext/src/gtest/gtest_main.cc -o CMakeFiles/gtest_main.dir/gtest_main.cc.s

CMakeFiles/gtest_main.dir/gtest_main.cc.o.requires:
.PHONY : CMakeFiles/gtest_main.dir/gtest_main.cc.o.requires

CMakeFiles/gtest_main.dir/gtest_main.cc.o.provides: CMakeFiles/gtest_main.dir/gtest_main.cc.o.requires
	$(MAKE) -f CMakeFiles/gtest_main.dir/build.make CMakeFiles/gtest_main.dir/gtest_main.cc.o.provides.build
.PHONY : CMakeFiles/gtest_main.dir/gtest_main.cc.o.provides

CMakeFiles/gtest_main.dir/gtest_main.cc.o.provides.build: CMakeFiles/gtest_main.dir/gtest_main.cc.o

# Object files for target gtest_main
gtest_main_OBJECTS = \
"CMakeFiles/gtest_main.dir/gtest_main.cc.o"

# External object files for target gtest_main
gtest_main_EXTERNAL_OBJECTS =

libgtest_main.a: CMakeFiles/gtest_main.dir/gtest_main.cc.o
libgtest_main.a: CMakeFiles/gtest_main.dir/build.make
libgtest_main.a: CMakeFiles/gtest_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libgtest_main.a"
	$(CMAKE_COMMAND) -P CMakeFiles/gtest_main.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gtest_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gtest_main.dir/build: libgtest_main.a
.PHONY : CMakeFiles/gtest_main.dir/build

CMakeFiles/gtest_main.dir/requires: CMakeFiles/gtest_main.dir/gtest_main.cc.o.requires
.PHONY : CMakeFiles/gtest_main.dir/requires

CMakeFiles/gtest_main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gtest_main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gtest_main.dir/clean

CMakeFiles/gtest_main.dir/depend:
	cd /ssd/ig_repertoire_constructor/ext/src/gtest && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ssd/ig_repertoire_constructor/ext/src/gtest /ssd/ig_repertoire_constructor/ext/src/gtest /ssd/ig_repertoire_constructor/ext/src/gtest /ssd/ig_repertoire_constructor/ext/src/gtest /ssd/ig_repertoire_constructor/ext/src/gtest/CMakeFiles/gtest_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gtest_main.dir/depend

