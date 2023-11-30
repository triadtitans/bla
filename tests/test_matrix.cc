#include <iostream>
#include "matrix.h"

namespace bla = ASC_bla;



int main()
{
  size_t n = 10;
  bla::Matrix<double> m(30,30);
  bla::Matrix<double> q(30,30);
  srand(7893549034785);

  for (int i = 0; i < 30; i++) {
    for (int j = 0; j < 30; j++) {
      m(i,j) = (int)(10*(rand()/(double)RAND_MAX));
      q(i,j) = (int)(10*(rand()/(double)RAND_MAX));
    }
  }
  std::cout << m << std::endl << std::endl;
  std::cout << q;

  //bla::MatrixView<double> view = m.Cols(0,2).Rows(1,1);
  //view = 1000;
  
  //m = bla::Matrix<double,bla::Ordering::RowMajor>(10.0 * m);





  //std::cout << "Inverse\n";
  //std::cout << q;
  //std::cout << q;
  //std::cout << bla::inverse(q);

  std::cout << std::endl << std::endl;
  std::cout << "fastMul" << std::endl;
  std::cout << m*q;
  std::cout << std::endl << std::endl;
  std::cout << fastMul(m,q);

  std::cout << std::endl << "difference" << std::endl;
  std::cout << (m*q) + (-1)*fastMul(m,q) << std::endl;
}