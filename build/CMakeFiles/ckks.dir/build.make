# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /opt/cmake-3.21.4/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.21.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /root/Homomorphic-Encryption-on-GPU

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /root/Homomorphic-Encryption-on-GPU/build

# Include any dependencies generated for this target.
include CMakeFiles/ckks.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ckks.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ckks.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ckks.dir/flags.make

CMakeFiles/ckks.dir/ckkstest.cu.o: CMakeFiles/ckks.dir/flags.make
CMakeFiles/ckks.dir/ckkstest.cu.o: ../ckkstest.cu
CMakeFiles/ckks.dir/ckkstest.cu.o: CMakeFiles/ckks.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Homomorphic-Encryption-on-GPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/ckks.dir/ckkstest.cu.o"
	/usr/local/cuda/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/ckks.dir/ckkstest.cu.o -MF CMakeFiles/ckks.dir/ckkstest.cu.o.d -x cu -dc /root/Homomorphic-Encryption-on-GPU/ckkstest.cu -o CMakeFiles/ckks.dir/ckkstest.cu.o

CMakeFiles/ckks.dir/ckkstest.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/ckks.dir/ckkstest.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/ckks.dir/ckkstest.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/ckks.dir/ckkstest.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target ckks
ckks_OBJECTS = \
"CMakeFiles/ckks.dir/ckkstest.cu.o"

# External object files for target ckks
ckks_EXTERNAL_OBJECTS =

CMakeFiles/ckks.dir/cmake_device_link.o: CMakeFiles/ckks.dir/ckkstest.cu.o
CMakeFiles/ckks.dir/cmake_device_link.o: CMakeFiles/ckks.dir/build.make
CMakeFiles/ckks.dir/cmake_device_link.o: CMakeFiles/ckks.dir/dlink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/Homomorphic-Encryption-on-GPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA device code CMakeFiles/ckks.dir/cmake_device_link.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ckks.dir/dlink.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ckks.dir/build: CMakeFiles/ckks.dir/cmake_device_link.o
.PHONY : CMakeFiles/ckks.dir/build

# Object files for target ckks
ckks_OBJECTS = \
"CMakeFiles/ckks.dir/ckkstest.cu.o"

# External object files for target ckks
ckks_EXTERNAL_OBJECTS =

ckks: CMakeFiles/ckks.dir/ckkstest.cu.o
ckks: CMakeFiles/ckks.dir/build.make
ckks: CMakeFiles/ckks.dir/cmake_device_link.o
ckks: CMakeFiles/ckks.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/Homomorphic-Encryption-on-GPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA executable ckks"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ckks.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ckks.dir/build: ckks
.PHONY : CMakeFiles/ckks.dir/build

CMakeFiles/ckks.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ckks.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ckks.dir/clean

CMakeFiles/ckks.dir/depend:
	cd /root/Homomorphic-Encryption-on-GPU/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/Homomorphic-Encryption-on-GPU /root/Homomorphic-Encryption-on-GPU /root/Homomorphic-Encryption-on-GPU/build /root/Homomorphic-Encryption-on-GPU/build /root/Homomorphic-Encryption-on-GPU/build/CMakeFiles/ckks.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ckks.dir/depend

