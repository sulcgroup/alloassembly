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
CMAKE_COMMAND = /snap/clion/203/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/203/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/josh/git/alloassembly

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/josh/git/alloassembly/build

# Utility rule file for evans.

# Include any custom commands dependencies for this target.
include contrib/evans/CMakeFiles/evans.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/evans/CMakeFiles/evans.dir/progress.make

contrib/evans/CMakeFiles/evans: ../contrib/evans/AllostericPatchySwapInteraction.so
contrib/evans/CMakeFiles/evans: ../contrib/evans/AllostericPatchyParticle.so
contrib/evans/CMakeFiles/evans: ../contrib/evans/AllostericPatch.so
contrib/evans/CMakeFiles/evans: ../contrib/evans/PLClusterTopology.so

evans: contrib/evans/CMakeFiles/evans
evans: contrib/evans/CMakeFiles/evans.dir/build.make
.PHONY : evans

# Rule to build all files generated by this target.
contrib/evans/CMakeFiles/evans.dir/build: evans
.PHONY : contrib/evans/CMakeFiles/evans.dir/build

contrib/evans/CMakeFiles/evans.dir/clean:
	cd /home/josh/git/alloassembly/build/contrib/evans && $(CMAKE_COMMAND) -P CMakeFiles/evans.dir/cmake_clean.cmake
.PHONY : contrib/evans/CMakeFiles/evans.dir/clean

contrib/evans/CMakeFiles/evans.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/contrib/evans /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/contrib/evans /home/josh/git/alloassembly/build/contrib/evans/CMakeFiles/evans.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/evans/CMakeFiles/evans.dir/depend
