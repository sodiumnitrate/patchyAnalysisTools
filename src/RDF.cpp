#define _USE_MATH_DEFINES
#include <cmath>
#include <pybind11/pybind11.h>
#include "pybind11/numpy.h"
#include <iostream>

namespace py = pybind11;

py::array_t<double> RDF_virtualcopies(
                                     py::array_t<double, py::array::c_style | py::array::forcecast> xcoords,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> ycoords,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> zcoords,
                                     double L,
                                     double cutoffvalue,
                                     double stepsize,
                                     const int nbins,
                                     double density,
                                     const int nparts)
{
    // process numpy array inputs from the python side
    auto buf2 = xcoords.request();
    double *x = (double *)buf2.ptr;

    auto buf3 = ycoords.request();
    double *y = (double *)buf3.ptr;

    auto buf4 = zcoords.request();
    double *z = (double *)buf4.ptr;

    // allocate mem for the output numpy array
    py::array_t<double> result = py::array_t<double>(nbins);
    auto buf5 = result.request();
    // c++ counterpart of the result array
    double *hist = (double *)buf5.ptr;

    // var declarations
    int ix, iy, iz, i, j, binnr;
    double dx, dy, dz, shellVol, r2;
    double c2 = cutoffvalue * cutoffvalue;


    // initialize result (hist) array to 0
    for (i = 0; i < nbins; i++)
    {
        hist[i] = 0;
    }

    // create histogram (the box with itself)
    for (i = 0; i < nparts; i++)
    {
        for (j = i+1; j < nparts; j++)
        {

            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];
            
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < c2)
            {
                binnr = int(sqrt(r2) / stepsize);
                hist[binnr] += 2;
            }
        }
    }

    // create histogram (the box with its neighbors)
    for (ix = -1; ix < 2; ix++)
    {
        for (iy = -1; iy < 2; iy++)
        {
            for (iz = -1; iz < 2; iz++)
            {
                if(ix == 0 && iy == 0 && iz == 0)
                    continue;
                for (i = 0; i < nparts; i++)
                {
                    for (j = 0; j < nparts; j++)
                    {
                        dx = x[j] + ix * L - x[i];
                        dy = y[j] + iy * L - y[i];
                        dz = z[j] + iz * L - z[i];

                        r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 < c2)
                        {
                            binnr = int(pow(r2, 0.5) / stepsize);
                            hist[binnr] += 1;
                        }
                    }
                }
            }
        }
    }

    //normalization
    for (i = 0; i < nbins; i++)
    {
        shellVol = ((pow(i + 1, 3.) - pow(i, 3.))) * pow(stepsize, 3.);
        hist[i] /= nparts * density * M_PI * 4. / 3 * shellVol;
    }

    return result;
}

py::array_t<double> RDF_virtualcopies_diff_grps(
                                     py::array_t<double, py::array::c_style | py::array::forcecast> xcoords_1,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> ycoords_1,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> zcoords_1,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> xcoords_2,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> ycoords_2,
                                     py::array_t<double, py::array::c_style | py::array::forcecast> zcoords_2,
                                     double L,
                                     double cutoffvalue,
                                     double stepsize,
                                     const int nbins,
                                     double density,
                                     const int nparts_tot,
                                     const int nparts1,
                                     const int nparts2)
{
    // process numpy array inputs from the python side
    auto buf2 = xcoords_1.request();
    double *x1 = (double *)buf2.ptr;

    auto buf3 = ycoords_1.request();
    double *y1 = (double *)buf3.ptr;

    auto buf4 = zcoords_1.request();
    double *z1 = (double *)buf4.ptr;

    auto buf5 = xcoords_2.request();
    double *x2 = (double *)buf5.ptr;

    auto buf6 = ycoords_2.request();
    double *y2 = (double *)buf6.ptr;

    auto buf7 = zcoords_2.request();
    double *z2 = (double *)buf7.ptr;

    // allocate mem for the output numpy array
    py::array_t<double> result = py::array_t<double>(nbins);
    auto buf8 = result.request();
    // c++ counterpart of the result array
    double *hist = (double *)buf8.ptr;

    // var declarations
    int ix, iy, iz, i, j, binnr;
    double dx, dy, dz, shellVol, r2;
    double c2 = cutoffvalue * cutoffvalue;

    // initialize result (hist) array to 0
    for (i = 0; i < nbins; i++)
    {
        hist[i] = 0;
    }

    // create histogram (the box with its neighbors)
    for (ix = -1; ix < 2; ix++)
    {
        for (iy = -1; iy < 2; iy++)
        {
            for (iz = -1; iz < 2; iz++)
            {
                for (i = 0; i < nparts1; i++)
                {
                    for (j = 0; j < nparts2; j++)
                    {
                        dx = x2[j] + ix * L - x1[i];
                        dy = y2[j] + iy * L - y1[i];
                        dz = z2[j] + iz * L - z1[i];

                        r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 < c2)
                        {
                            binnr = int(pow(r2, 0.5) / stepsize);
                            hist[binnr] += 1;
                        }
                    }
                }
            }
        }
    }

    //normalization
    for (i = 0; i < nbins; i++)
    {
        shellVol = ((pow(i + 1, 3.) - pow(i, 3.))) * pow(stepsize, 3.);
        //hist[i] /= (nparts_tot/((float)nparts1*nparts2)) * density * M_PI * 4. / 3 * shellVol;
        hist[i] /= (nparts1*nparts2/(float)nparts_tot) * density * M_PI * 4. / 3 * shellVol;
    }

    return result;
}

PYBIND11_MODULE(RDF, m)
{
    m.doc() = "two methods for calculating radial distribution functions";
    m.def("RDF_virtualcopies", &RDF_virtualcopies, "calculate RDF with virtual copies");
    m.def("RDF_virtualcopies_diff_grps", &RDF_virtualcopies_diff_grps, "calculate RDF with virtual copies");
}
