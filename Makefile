# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/DATA/Document/code_exp/cpp_example/MatrixMultiply/c_array

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/DATA/Document/code_exp/cpp_example/MatrixMultiply/c_array

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /media/DATA/Document/code_exp/cpp_example/MatrixMultiply/c_array/CMakeFiles /media/DATA/Document/code_exp/cpp_example/MatrixMultiply/c_array/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /media/DATA/Document/code_exp/cpp_example/MatrixMultiply/c_array/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named matspeedtest

# Build rule for target.
matspeedtest: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 matspeedtest
.PHONY : matspeedtest

# fast build rule for target.
matspeedtest/fast:
	$(MAKE) -f CMakeFiles/matspeedtest.dir/build.make CMakeFiles/matspeedtest.dir/build
.PHONY : matspeedtest/fast

#=============================================================================
# Target rules for targets named mattest

# Build rule for target.
mattest: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mattest
.PHONY : mattest

# fast build rule for target.
mattest/fast:
	$(MAKE) -f CMakeFiles/mattest.dir/build.make CMakeFiles/mattest.dir/build
.PHONY : mattest/fast

src/test/matspeedtest.o: src/test/matspeedtest.cpp.o
.PHONY : src/test/matspeedtest.o

# target to build an object file
src/test/matspeedtest.cpp.o:
	$(MAKE) -f CMakeFiles/matspeedtest.dir/build.make CMakeFiles/matspeedtest.dir/src/test/matspeedtest.cpp.o
.PHONY : src/test/matspeedtest.cpp.o

src/test/matspeedtest.i: src/test/matspeedtest.cpp.i
.PHONY : src/test/matspeedtest.i

# target to preprocess a source file
src/test/matspeedtest.cpp.i:
	$(MAKE) -f CMakeFiles/matspeedtest.dir/build.make CMakeFiles/matspeedtest.dir/src/test/matspeedtest.cpp.i
.PHONY : src/test/matspeedtest.cpp.i

src/test/matspeedtest.s: src/test/matspeedtest.cpp.s
.PHONY : src/test/matspeedtest.s

# target to generate assembly for a file
src/test/matspeedtest.cpp.s:
	$(MAKE) -f CMakeFiles/matspeedtest.dir/build.make CMakeFiles/matspeedtest.dir/src/test/matspeedtest.cpp.s
.PHONY : src/test/matspeedtest.cpp.s

src/test/mattest.o: src/test/mattest.cpp.o
.PHONY : src/test/mattest.o

# target to build an object file
src/test/mattest.cpp.o:
	$(MAKE) -f CMakeFiles/mattest.dir/build.make CMakeFiles/mattest.dir/src/test/mattest.cpp.o
.PHONY : src/test/mattest.cpp.o

src/test/mattest.i: src/test/mattest.cpp.i
.PHONY : src/test/mattest.i

# target to preprocess a source file
src/test/mattest.cpp.i:
	$(MAKE) -f CMakeFiles/mattest.dir/build.make CMakeFiles/mattest.dir/src/test/mattest.cpp.i
.PHONY : src/test/mattest.cpp.i

src/test/mattest.s: src/test/mattest.cpp.s
.PHONY : src/test/mattest.s

# target to generate assembly for a file
src/test/mattest.cpp.s:
	$(MAKE) -f CMakeFiles/mattest.dir/build.make CMakeFiles/mattest.dir/src/test/mattest.cpp.s
.PHONY : src/test/mattest.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... matspeedtest"
	@echo "... mattest"
	@echo "... rebuild_cache"
	@echo "... src/test/matspeedtest.o"
	@echo "... src/test/matspeedtest.i"
	@echo "... src/test/matspeedtest.s"
	@echo "... src/test/mattest.o"
	@echo "... src/test/mattest.i"
	@echo "... src/test/mattest.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
