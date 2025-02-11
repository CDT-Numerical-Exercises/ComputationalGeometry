#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"
#include "convex_hull.h"
#include "csv.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Requires path to CSV file as arg." << std::endl;
    return 1;
  }

  gsl_matrix *data = load_csv_to_dmatrix(argv[1]);

  std::vector<size_t> hull;
  if (argc < 3) {
    hull = convex_hull(data);
  } else {
    hull = convex_hull(data, argv[2]);
  }

  std::cout << "Hull contains " << hull.size() << " points" << std::endl;

  Gnuplot gp;
  gp << "set key off\n";
  // gp << "set yrange[-0.5:2]\nset xrange[-0.5:1.5]\n";
  gp << "plot '-' with circles, "; // for the points
  gp << "'-' with lines\n"; // for the hull

  // plot the points
  for (int row = 0; row < data->size1; ++row) {
    for (int col = 0; col < data->size2 - 1; ++col) {
      gp << gsl_matrix_get(data, row, col) << " ";
    }
    gp << gsl_matrix_get(data, row, data->size2 - 1);
    gp << "\n";
  }
  gp << "e\n";

  // draw the hull
  for (const size_t &i : hull) {
    gsl_vector_const_view row_view = gsl_matrix_const_row(data, i);
    gp << gsl_vector_get(&row_view.vector, 0) << " " << gsl_vector_get(&row_view.vector, 1) << "\n";
  }
  {
    gsl_vector_const_view row_view = gsl_matrix_const_row(data, hull[0]);
    gp << gsl_vector_get(&row_view.vector, 0) << " " << gsl_vector_get(&row_view.vector, 1) << "\n";
  }
  gp << "e\n";

  gsl_matrix_free(data);
}
