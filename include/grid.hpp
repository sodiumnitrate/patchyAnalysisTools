/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include <vector>
#include <Eigen/Dense>
#include <iostream>

template <class T>
class Grid{
    int x;
    int y;
    int z;
    std::vector<std::vector<std::vector<T> > > values;
public:
    // constructor
    Grid(int x_, int y_, int z_, T val);

    // getters
    std::vector<int> get_size();
    T operator()(int ix, int iy, int iz);
    void print();

    // sum operations
    std::vector<std::vector<T> > sum(int axis);
    std::vector<std::vector<T> > sum_along_x();
    std::vector<std::vector<T> > sum_along_y();
    std::vector<std::vector<T> > sum_along_z();
    void add_subgrid(Grid<T> subgrid, int lx, int ly, int lz);
};