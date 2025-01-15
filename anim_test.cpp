#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"

int main() {
  Gnuplot gp;

  std::filesystem::create_directory("output");

  gp << "set terminal pngcairo size 350,262 enhanced font 'Verdana,10'\n";
  gp << "set output 'output/frame001.png'\n";
  gp << "set yrange[-1:1]\n";
  gp << "plot '-' with circles\n";
  gp << "0.0 0.0\ne\n";
  gp << "set terminal pngcairo size 350,262 enhanced font 'Verdana,10'\n";
  gp << "set output 'output/frame002.png'\n";
  gp << "plot '-' with circles\n";
  gp << "0.0 0.0\n";
  gp << "0.0 1.0\n";

  return 0;
}
