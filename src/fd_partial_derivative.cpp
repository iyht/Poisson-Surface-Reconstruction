#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;
    tripleList.reserve(nx*ny*nz);

    int new_nx = nx, new_ny = ny, new_nz = nz;
    int x_dir = 0, y_dir = 0, z_dir = 0;

    if (dir == 0)
    {
        new_nx--; x_dir = 1;
    }
    else if (dir == 1)
    {
        new_ny--; y_dir = 1;
    }
    else
    {
        new_nz--; z_dir = 1;
    }

    // loop through 3 dims
    for (int i = 0; i < new_nx; i++) { // x
        for (int j = 0; j < new_ny; j++) { // y
            for (int k = 0; k < new_nz; k++) { //z
                // Note: D: #eqs by nx*ny*nz
                tripleList.push_back(T(i + new_nx * (j + k * new_ny), i + nx * (j + k * ny), -1.0 / h));
                tripleList.push_back(T(i + new_nx * (j + k * new_ny), i + x_dir + nx * ((j + y_dir) + (k + z_dir) * ny), 1.0 / h));
            }
        }
    }
    D.resize(new_nx*new_ny*new_nz, nx*ny*nz);
    D.setFromTriplets(tripleList.begin(), tripleList.end());

}
