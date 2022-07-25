#define _USE_MATH_DEFINES
#include <cmath>
#include <pybind11/pybind11.h>
#include "pybind11/numpy.h"
#include <iostream>

namespace py = pybind11;

py::array_t<double> sq(const int g, double L, py::array_t<double, py::array::c_style | py::array::forcecast> points, int numpoints)
{
    // process numpy array inputs from the python side
    auto buf1 = points.request();
    double *pts = (double *)buf1.ptr;

    // allocate mem for the output numpy array
    py::array_t<double> result = py::array_t<double>(2 * g);
    auto buf3 = result.request();

    // c++ counterpart of the numpy array output
    double *sqPtr = (double *)buf3.ptr;

    // var declarations
    int i, ix, iy, iz;
    double qx;
    double qy;
    double qz;
    double *c = new double[2 * g];
    double *s = new double[2 * g];
    int *count = new int[2 * g];
    double binsize = 2 * M_PI / L;
    double qval, cval, sval, ptx, pty, ptz, sSum, cSum, rdotq, npts, ct;
    int qind;

    // initialize vars
    for (i = 0; i < 2 * g; i++)
    {
        c[i] = 0;
        s[i] = 0;
        count[i] = 0;
        sqPtr[i] = 0;
    }

    // loop over all grid points
    for (ix = -g; ix <= g; ix++)
    {
        qx = ix * binsize;
        for (iy = -g; iy <= g; iy++)
        {
            qy = iy * binsize;
            for (iz = -g; iz <= g; iz++)
            {
                qz = iz * binsize;

                // norm of q = (qx, qy, qz)
                qval = pow(qx * qx + qy * qy + qz * qz, 0.5);

                // bin index for the given q vector
                qind = (int)(qval / (double)binsize);

                // increment count for the given index
                count[qind] += 1;

                cSum = 0.;
                sSum = 0.;
                for (i = 0; i < numpoints; i++)
                {
                    // get point coords
                    ptx = pts[i];
                    pty = pts[numpoints + i];
                    ptz = pts[2*numpoints + i];

                    // r dot q
                    rdotq = ptx * qx + pty * qy + ptz * qz;

                    // cos(r dot q)
                    cval = cos(rdotq);

                    // sin(r dot q)
                    sval = sin(rdotq);

                    // sum cosine values for the given index
                    cSum += cval;

                    // sum sine values for the given index
                    sSum += sval;
                }

                // square sum of cosines
                //cSum = cSum * cSum;
                cSum = pow(cSum,2.);

                // square sum of sines
                //sSum = sSum * sSum;
                sSum = pow(sSum,2.);

                c[qind] += cSum;
                s[qind] += sSum;
            }
        }
    }

    npts = double(numpoints);

    // calculate sq
    for (i = 0; i < 2 * g; i++)
    {
        ct = double(count[i]);
        if (count[i] > 0)
            sqPtr[i] = (1./ npts) * ((c[i] + s[i]) / ct);
        else
            sqPtr[i] = 0.;
    }

    return result;
}

PYBIND11_MODULE(sq, m)
{
    m.doc() = "calculates s(q) from coordinates";
    m.def("sq", &sq, "A function that calculates s(q) given coordinates");
}
