# AutoSync

AutoSync is a framework which provides automatic generation of synchronisation primitives for multi-threaded C programs. The framework consists in a simple descriptive language that enables developers to express their synchronisation intentions, instead of coding synchronisation explicitly. AutoSync is then able to statically analyse these intentions to generate synchronisation mechanisms that prevent concurrency bugs, while assuring optimal performance. Moreover, the tool facilitates the understanding of synchronization among threads, provides necessary information for designing test cases, and make the concurrency-specific parts of source code more readable. Based on the explicit synchronisation intentions, it can be used to enable visualization of concurrent software properties.

![AutoSync Overview](doc/AutoSyncOverview.png)

# How it works
First, developers write their code using the provided descriptive language (see [AutoSync.h](src/AutoSync.h)). Then, the framework is able to parse the original code, by means of Static Analysis, and extract the highest number of program features that can be later on used to fine tune the generated code. The Code Generator itself uses this information with some developed algorithms in order to write the new code and replace it in the original files before compilation. Thus, the Code Generator delivers new files that contain only synchronization primitives from the standard libraries and are ready to be compiled. The compilation itself is not part of the framework. AutoSync receive source code files and produces as output new source files. That way, the process of compilation and linking is left open to be performed as desired.

# How to use it
The following steps are required when writing code using the AutoSync framework:

1. Include the file "AutoSync.h" in our project;
2. Express the desired high-level synchronisation intentions using the methods present in the [interface](src/AutoSync.h);
3. Call the [Makefile](Makefile) providing the path to your C file;
   `````
   $ make auto_autosync some_file.c 
   `````
4. In the folder [generated](generated/) you can find the new files with the generated code, which can be built as usual.
[WIP]

# Benchmarks
The FFT program from the well-known SPLASH benchmark has been refactored to evaluate AutoSync. The original version can be found [here](https://github.com/SakalisC/Splash-3/blob/master/codes/kernels/fft/fft.c.in).
The refactored version is [here](examples/benchmark_splash_fft/fft_auto_sync.c).


# Contributors
AutoSync was developed during the Master Thesis program of Software Engineering for Embedded Systems in the Rheinland-Pfälzische Technische Universität Kaiserslautern-Landau (Germany) by Matheus Bortoloti under the supervision of Dr. Jasmin Jahić.
