# CircuitSimulator_MAI_JUK1

This repo showcases the simulation code developed as part of my MAI final project. It is an RF simulator with extensibility built in mind.

## Requirements

This project requires a compiler capable of compiling C++20 or higher

For vector fitting functionality Matlab is required to perform the fit. Otherwise a PRR file is need ([See PRR files for an example](src/PRR/PRRFile Format.md))

For outputting graphs directly Python3 support with matplotlib and numpy is required

## Compiling

This project has only been compiled on windows, but should theoretically compile on other operating systems. clang 10.0.0 was used for testing. The pop-up window for file selection is only available on windows. A bat file for windows has been provided for building in both release and debug mode

## Netlist Files

Information on netlist files is provided in the netlist folder. Note the units of the simulation are non-standard

## Running

The simulator can be launched on windows without a target, it will prompt the user to select a netlist file. On other operating systems, one must be sure pass in the desired netlist to the program, and ensure the working directory is compatible with all paths used in the netlist

## Documentation

Code documentation is provided via the included PDF [DoxygenDocumentation.pdf](DoxygenDocumentation.pdf). A website version is also hosted [on my github pages here](https://dslogget.github.io/JUK1/index.html)

The final produced MAI thesis is provided in the PDF [DSloggettSchEng2021.pdf](DSloggettSchEng2021.pdf)

## Acknowledgements

I would like to thank my master's supervisor Dr. Justin King for helping me grasp the concepts involved in this simulation throughout the project. Without his guidance, I would not have been able to come close to understanding how RF simulations work
