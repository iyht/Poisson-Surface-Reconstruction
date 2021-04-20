#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here

  int m = P.rows();
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripleList;
  tripleList.reserve(8*m);
  for(int r = 0; r < m; r++)
  {
      // initial the value for (x, y, z) of a point and get the position w.r.t the corner
      double p_x = P(r, 0) - corner(0);
      double p_y = P(r, 1) - corner(1);
      double p_z = P(r, 2) - corner(2);

      // get position of the corner of the cube where the point is located w.r.t grid.
      int i = (int)floor(p_x/h);
      int j = (int)floor(p_y/h);
      int k = (int)floor(p_z/h);

      // get the position of point w.r.t to the cube.
      double c_x = (p_x - i*h)/h;
      double c_y = (p_y - j*h)/h;
      double c_z = (p_z - k*h)/h;

      // calculate the trilinear interpolation weight
      tripleList.push_back(T(r, i+(nx*j)+(ny*nx*k), (1-c_x)*(1-c_y)*(1-c_z)));
      tripleList.push_back(T(r, (i+1)+(nx*j)+(ny*nx*k), (c_x)*(1-c_y)*(1-c_z)));
      tripleList.push_back(T(r, i+(nx*(j+1))+(ny*nx*k), (1-c_x)*(c_y)*(1-c_z)));
      tripleList.push_back(T(r, i+(nx*j)+(ny*nx*(k+1)), (1-c_x)*(1-c_y)*(c_z)));
      tripleList.push_back(T(r, (i+1)+(nx*j)+(ny*nx*(k+1)), (c_x)*(1-c_y)*(c_z)));
      tripleList.push_back(T(r, i+(nx*(j+1))+(ny*nx*(k+1)), (1-c_x)*(c_y)*(c_z)));
      tripleList.push_back(T(r, (i+1)+(nx*(j+1))+(ny*nx*k), (c_x)*(c_y)*(1-c_z)));
      tripleList.push_back(T(r, (i+1)+(nx*(j+1))+(ny*nx*(k+1)), (c_x)*(c_y)*(c_z)));

  }
  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(tripleList.begin(), tripleList.end());


    ////////////////////////////////////////////////////////////////////////////
}
