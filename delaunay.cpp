#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <unordered_set>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"
#include "delaunay.h"

const double data[4*2] = {
  0, 0,
  1, 1,
  -0.25, 1.6,
  1, -0.25
};

int main(int argc, char *argv[]) {
  const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 4, 2);

  std::vector<Triangle> triangulation;
  if (argc < 2) {
    triangulation = delaunay_triangulate(&d_view.matrix);
  } else {
    triangulation = delaunay_triangulate(&d_view.matrix, argv[1]);
  }

  std::cout << "Triangulation produced " << triangulation.size() << " triangles." << std::endl;

  // add the edges to an unordered set
  // this will make plotting easier
  std::unordered_set<Edge> edges;
  for (const Triangle &tri : triangulation) {
    Edge A; Edge B; Edge C;
    tri.to_edges(A, B, C);
    edges.insert(A);
    edges.insert(B);
    edges.insert(C);
  }

  // set up gnuplot
  Gnuplot gp;
  gp << "set key off\n";
  gp << "set yrange[-0.5:2]\nset xrange[-0.5:1.5]\n";
  gp << "plot '-' with circles, "; // for the points

  // add a curve for every edge
  const size_t N_edges = edges.size();
  for (size_t i = 1; i < N_edges; ++i) {
    gp << "'-' with lines, ";
  }
  gp << "'-' with lines\n";

  // plot the points
  for (int row = 0; row < d_view.matrix.size1; ++row) {
    for (int col = 0; col < d_view.matrix.size2 - 1; ++col) {
      gp << gsl_matrix_get(&d_view.matrix, row, col) << " ";
    }
    gp << gsl_matrix_get(&d_view.matrix, row, d_view.matrix.size2 - 1);
    gp << "\n";
  }
  gp << "e\n";

  // draw the lines
  for (const Edge &e : edges) {
    gsl_vector_const_view A_view = gsl_matrix_const_row(&d_view.matrix, e.A);
    gp << gsl_vector_get(&A_view.vector, 0) << " " << gsl_vector_get(&A_view.vector, 1) << "\n";
    gsl_vector_const_view B_view = gsl_matrix_const_row(&d_view.matrix, e.B);
    gp << gsl_vector_get(&B_view.vector, 0) << " " << gsl_vector_get(&B_view.vector, 1) << "\ne\n";
  }
}
