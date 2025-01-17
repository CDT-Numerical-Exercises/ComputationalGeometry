#include <iostream>
#include <unordered_set>

#include "delaunay.h"

int main() {
  const Edge A(1, 2);
  const Edge B(2, 1);
  const Edge C(3, 1);

  std::unordered_set<Edge> s;
  s.insert(A);
  s.insert(B);
  s.insert(C);

  for (Edge e : s) {
    std::cout << "Edge(" << e.A << ", " << e.B << ")" << std::endl;
  }

  return 0;
}
