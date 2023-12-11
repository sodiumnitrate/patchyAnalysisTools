/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#include "grid.hpp"

// constructor
template <class T>
Grid<T>::Grid(int x_, int y_, int z_, T val){
    x = x_;
    y = y_;
    z = z_;
    for(int i = 0; i < x; i++){
        std::vector<std::vector<T> > plane;
        for(int j = 0; j < y; j++){
            std::vector<T> arr;
            for(int k = 0; k < z; k++){
                arr.push_back(val);
            }
            plane.push_back(arr);
        }
        values.push_back(plane);
    }
}

// getters
template <class T> std::vector<int> Grid<T>::get_size(){
    std::vector<int> size;
    size.push_back(x);
    size.push_back(y);
    size.push_back(z);
    return size;
}
template <class T> T Grid<T>::operator()(int ix, int iy, int iz){
    return values[ix][iy][iz];
}
template <class T> void Grid<T>::print(){
    for(int i = 0; i < z; i++){
        std::cout << "z = " << i << std::endl;
        for(int j = 0; j < x; j++){
            for(int k = 0; k < y; k++){
                std::cout << values[j][k][i] << ", ";
            }
            std::cout << std::endl;
        }
    }
}

// reduce dim
template <class T> std::vector<std::vector<T> > Grid<T>::sum(int axis){
    switch(axis){
        case 0:
            return this->sum_along_x();
            break;
        case 1:
            return this->sum_along_y();
            break;
        case 2:
            return this->sum_along_z();
            break;
        default:
            throw "invalid axis value";
    }
}

template <class T> std::vector<std::vector<T> > Grid<T>::sum_along_x(){
    std::vector<std::vector<T> > plane;
    plane.resize(y, std::vector<T>(z));
    for (auto &i : plane) std::fill(i.begin(), i.end(), 0);
    for (int iy = 0; iy < y; iy++){
        for(int iz = 0; iz < z; iz++){
            for(int ix = 0; ix < x; ix++){
                plane[iy][iz] += values[ix][iy][iz];
            }
        }
    }
    return plane;
}

template <class T> std::vector<std::vector<T> > Grid<T>::sum_along_y(){
    std::vector<std::vector<T> > plane;
    plane.resize(x, std::vector<T>(z));
    for (auto &i : plane) std::fill(i.begin(), i.end(), 0);
    for (int ix = 0; ix < x; ix++){
        for(int iz = 0; iz < z; iz++){
            for(int iy = 0; iy < y; iy++){
                plane[ix][iz] += values[ix][iy][iz];
            }
        }
    }
    return plane;
}

template <class T> std::vector<std::vector<T> > Grid<T>::sum_along_z(){
    std::vector<std::vector<T> > plane;
    plane.resize(x, std::vector<T>(y));
    for (auto &i : plane) std::fill(i.begin(), i.end(), 0);
    for (int ix = 0; ix < x; ix++){
        for(int iy = 0; iy < y; iy++){
            for(int iz = 0; iz < z; iz++){
                plane[ix][iy] += values[ix][iy][iz];
            }
        }
    }
    return plane;
}

// add subgrid
template <class T> void Grid<T>::add_subgrid(Grid<T> subgrid, int lx, int ly, int lz){
    /*
    Adds elements of the subgrid at (lx,ly,lz). Wraps indices if need be.
    */
   int cx, cy, cz;
   std::vector<int> size_subgrid = subgrid.get_size();
   for(int ix = 0; ix < size_subgrid[0]; ix++){
        for(int iy = 0; iy < size_subgrid[1]; iy++){
            for(int iz = 0; iz < size_subgrid[2]; iz++){
                cx = ix + lx;
                cy = iy + ly;
                cz = iz + lz;
                if ( cx >= x ) cx -= x;
                else if ( cx < 0 ) cx += x;
                if ( cy >= y ) cy -= y;
                else if ( cy < 0 ) cy += y;
                if ( cz >= z ) cz -= z;
                else if ( cz < 0 ) cz += z;

                values[cx][cy][cz] += subgrid(ix,iy,iz);
            }
        }
   }
}

// allowed types
template class Grid<int>;
template class Grid<float>;
template class Grid<double>;