# HPC Project
![Build Status](https://img.shields.io/badge/build-passing-brightgreen)

## Building the Project

### Release Mode

To build the project in release mode, run the following commands:

```sh
./build.sh
```

This script will create a `build` directory, run CMake with the `Release` build type, and build the project.

## Floating Point Precision

The precision could be set with `-DPRECISION=float`, or simply use the shell script

```sh
./build_float.sh
```

This script will create a `build_float` directory, run CMake with the `Release` build type, and build the project.


### Debug Mode

To build the project in debug mode, run the following commands:

```sh
./build_debug.sh
```

This script will create a `build_debug` directory, run CMake with the `Debug` build type, and build the project.

## Running Tests

To run the tests, follow these steps:

1. Build the project (either in release or debug mode).
2. Navigate to the `build/tests` directory.
3. Run the test executable

```sh
./tests
```
or with ctest:
```sh
ctest --output-on-failure
```

The tests will be executed, and the results will be displayed in the terminal.

## Parameter File

The parameter file examples are stored in `examples/`, it should be self-explanatory.

The string function is parsed with C++ Mathematical Expression Toolkit Library [exprtk](https://github.com/ArashPartow/exprtk/tree/master).

## Dependencies

This project uses the following dependencies:

- [Catch2](https://github.com/catchorg/Catch2): A unit testing framework for C++.
- [exprtk](https://github.com/ArashPartow/exprtk/tree/master): C++ Mathematical Expression Toolkit Library

The dependencies are automatically fetched and built using CMake.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.