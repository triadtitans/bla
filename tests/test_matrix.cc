#include <iostream>

#include "matrix.h"

namespace bla = ASC_bla;

int main()
{
  size_t n = 10;
  bla::Matrix<int> m(3,3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m(i,j) = i+j;
    }
  }

  bla::MatrixView<int> view = m.Cols(0,2).Rows(1,1);
  
  view = 1000;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%d ", m(i,j));
    }
    printf("\n");
  }
}