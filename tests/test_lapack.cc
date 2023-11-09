#include <iostream>

#include <vector.h>
#include <lapack_interface.h>
#include <chrono>
#include <cmath>
#include <stdlib.h>

using namespace ASC_bla;
using namespace std;


int main()
{
  Vector<double> x(5);
  Vector<double> y(5);

Matrix<double> q(3,3);
Matrix<double> p(3,3);
Matrix<double> c(3,3);
  q(0,0)=1;
  q(0,1)=2;
  q(0,2)=3;
  q(1,1)=4;
  q(1,2)=5;
  q(2,2)=6;
  p(0,0)=1;
  p(0,1)=2;
  p(0,2)=3;
  p(1,1)=4;
  p(1,2)=5;
  p(2,2)=6;
  for (int i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 2;
    }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  
  AddVectorLapack (2, x, y);  
  cout << "y+2*x = " << y << endl;

  MultMatMatLapack(p,q,c);
  std::cout << c;
  srand(0);
  for(int n=1; n<=1;n++){
    Matrix<double> m1(pow(10,n),pow(10,n));
    Matrix<double> m2(pow(10,n),pow(10,n));
    Matrix<double> m3(pow(10,n),pow(10,n));
    for(int i=0;i<pow(10,n);i++){
      for(int j=0;j<pow(10,n);j++){
        m1(i,j)=(double)(rand())/RAND_MAX;
        m2(i,j)=(double)(rand())/RAND_MAX;
        m3(i,j)=0;
      }
    }
    size_t flops = n*n*n;
    size_t runs = size_t (1e9 / flops) + 1;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < runs; i++){
      MultMatMatLapack(m1,m2,m3);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end-start).count();
          
    cout << "n = " << n << ", time = " << time << " s, GFlops = " 
        << (n*runs)/time*1e-9 << endl;
  }

}

  
