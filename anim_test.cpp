#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"

const std::string dirname("animtest_output");

void save_plot(const std::filesystem::path fn, const std::string gpscript,
               const gsl_matrix *data) {
  Gnuplot gp;
  gp << "set terminal pngcairo size 350,262 enhanced font 'Verdana,10'\n";
  gp << "set output " << fn << "\n"; // be warned -- this may not be sanitised
  gp << gpscript << "\n";
  for (int row = 0; row < data->size1; ++row) {
    for (int col = 0; col < data->size2 - 1; ++col) {
      gp << gsl_matrix_get(data, row, col) << " ";
    }
    gp << gsl_matrix_get(data, row, data->size2 - 1);
    gp << "\n";
  }
  gp << "e\n";
}

int main() {
  Gnuplot gp;

  std::filesystem::create_directory(dirname);

  const double points[] =
  { 0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0,
    -1.0, 0.0,
    0.0, -1.0 };

  const gsl_matrix_const_view all_data = gsl_matrix_const_view_array(points, 5, 2);

  for (int points = 1; points <= 5; ++points) {
    const gsl_matrix_const_view slice_view = gsl_matrix_const_submatrix(&all_data.matrix, 0, 0, points, 2);
    char buf[50];
    snprintf(buf, sizeof(buf), (dirname+"/frame%03d.png").c_str(), points);
    std::cout << buf << std::endl;
    save_plot(buf, "set yrange[-1.5:1.5]\nset xrange[-1.5:1.5]\nplot '-' with circles\n", &slice_view.matrix);
  }

  // frames can then be merged together into an animation using ffmpeg
  // see: https://trac.ffmpeg.org/wiki/Slideshow
  /* something like
     ffmpeg -framerate 1 -i frame%03d.png -vf format=yuv420p -vcodec libx264 -preset veryslow -tune stillimage anim.mp4
     should do the trick */

  return 0;
}
