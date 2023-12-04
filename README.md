# patchyAnalysisTools

(In the process of rewriting the whole thing in C++, with python bindings, to support efficient calculation of various trajectory-wide quantities).

## TODO:
1. GSD C API. Handling things on the python side for now.
2. Proper error handling from the C++ side, visible to python as well. Using `throw`s for now.
3. ~~C++-level tests. Sticking to `pytest` for now.~~
4. Libraries for rotations.