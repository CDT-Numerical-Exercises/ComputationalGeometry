#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "convex_hull.h"

const double data[10*2] = {
  0.2, 0.85,
  0.3, 0.83,
  0.32, 0.79,
  0.42, 0.62,
  0.02, 0.52,
  0.415, 0.52,
  0.62, 0.52,
  0.41, 0.2,
  0.55, 0.12,
  0.27, 0.09
};

int main(int argc, char *argv[]) {
  const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 10, 2);

  std::vector<size_t> hull;
  if (argc < 2) {
    hull = convex_hull(&d_view.matrix);
  } else {
    hull = convex_hull(&d_view.matrix, argv[1]);
  }

  std::cout << "Hull contains " << hull.size() << " points" << std::endl;

  Gnuplot gp;
  gp << "set key off\n";
  // gp << "set yrange[-0.5:2]\nset xrange[-0.5:1.5]\n";
  gp << "plot '-' with circles, "; // for the points
  gp << "'-' with lines\n"; // for the hull

  // plot the points
  for (int row = 0; row < d_view.matrix.size1; ++row) {
    for (int col = 0; col < d_view.matrix.size2 - 1; ++col) {
      gp << gsl_matrix_get(&d_view.matrix, row, col) << " ";
    }
    gp << gsl_matrix_get(&d_view.matrix, row, d_view.matrix.size2 - 1);
    gp << "\n";
  }
  gp << "e\n";

  // draw the hull
  for (const size_t &i : hull) {
    gsl_vector_const_view row_view = gsl_matrix_const_row(&d_view.matrix, i);
    gp << gsl_vector_get(&row_view.vector, 0) << " " << gsl_vector_get(&row_view.vector, 1) << "\n";
  }
  {
    gsl_vector_const_view row_view = gsl_matrix_const_row(&d_view.matrix, hull[0]);
    gp << gsl_vector_get(&row_view.vector, 0) << " " << gsl_vector_get(&row_view.vector, 1) << "\n";
  }
  gp << "e\n";
}
