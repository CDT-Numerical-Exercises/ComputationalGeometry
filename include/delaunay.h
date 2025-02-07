#ifndef DELAUNAY_H
#define DELAUNAY_H 1

#include <vector>
#include <filesystem>
#include <gsl/gsl_matrix.h>

struct Edge {
  Edge() : A(0), B(0) {}
  Edge(const Edge &e) : A(e.A), B(e.B) {}
  Edge(const size_t a, const size_t b) : A(a), B(b) {}
  ~Edge() {}
  bool operator==(const Edge b) const;
  size_t A;
  size_t B;
};

// Specialize std::hash for Edge
// this way, we can use it without needing to include boost there
std::size_t hash_value(Edge const &p);
namespace std {
  template <>
  struct hash<Edge> {
    size_t operator()(const Edge& edge) const {
      return hash_value(edge);
    }
  };
}

class Triangle {
public:
  Triangle(const size_t a, const size_t b, const size_t c) : A(a), B(b), C(c) {}
  ~Triangle() {}

  int to_matrix(const gsl_matrix *verts, const gsl_matrix *supertriangle, gsl_matrix *X) const;
  int to_matrix(const gsl_matrix *verts, gsl_matrix *X) const;
  void to_edges(Edge &a, Edge &b, Edge &c) const;
  
  size_t A;
  size_t B;
  size_t C;

private:
  bool is_in_circumcircle(const gsl_vector *point, const gsl_matrix *verts, const gsl_matrix *supertriangle) const;

  friend std::vector<Triangle> delaunay_triangulate(const gsl_matrix *verts, const std::filesystem::path frame_dir, const bool save_anim);
};

int make_supertriangle(const gsl_matrix *data, gsl_matrix *supertriangle, const double dscale = 100);
int find_circumcircle(const gsl_matrix *points, double &r, gsl_vector *X);

std::vector<Triangle> delaunay_triangulate(const gsl_matrix *verts);
std::vector<Triangle> delaunay_triangulate(const gsl_matrix *verts, const std::filesystem::path frame_dir);

#endif
