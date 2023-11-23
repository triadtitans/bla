#include <iostream>
#include "matrix.h"

namespace bla = ASC_bla;



int main()
{
  size_t n = 10;
  bla::Matrix<double> m(3,3);
  bla::Matrix<double> q(3,3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m(i,j) = i+j;
    }
  }
  q(0,0)=1;
  q(0,1)=2;
  q(0,2)=3;
  q(1,1)=4;
  q(1,2)=5;
  q(2,2)=6;
  std::cout << m;
  bla::MatrixView<double> view = m.Cols(0,2).Rows(1,1);
  
  view = 1000;
  
  m = bla::Matrix<double,bla::Ordering::RowMajor>(10.0 * m);


  std::cout << q;


  std::cout << "Inverse\n";
  std::cout << q;
  std::cout << q;
  std::cout << bla::inverse(q);

  std::cout << m*q;
  std::cout << fastMul(m,q);
}