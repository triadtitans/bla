#include <iostream>

#include <vector.h>
#include "matrix.h"

namespace bla = ASC_bla;


int main()
{
  size_t n = 5;
  bla::Vector<double> x(n), y(n);
  bla::Matrix<double> a(n,n), b(n,n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  a.Dump();

  bla::Vector<double> z = x+y;
  
  std::cout << "x+y = " << z << std::endl;
}


