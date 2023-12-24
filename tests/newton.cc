#include <iostream>
#include "matrix.h"
#include "vector.h"
#include <functional>

using namespace ASC_bla;

class FiniteDifference{
    int N=0;
    double h=0;
    std::function<double(double)> f;
public:
    FiniteDifference(std::function<double(double)> func,int steps):N{steps},f{func}{
        if(steps < 2)
            throw std::invalid_argument("At least two steps required");
        h = 1.0/(1+N);
    }
    Vector<double> Eval(Vector<double> y){
        Vector<double> result(N);
        result(0)=2*y(0)-y(1)+h*h*(pow(y(0),3)-f(X(0)));
        result(N-1)=-y(N-2)+2*y(N-1)+h*h*(pow(y(N-1),3)-f(X(N-1)));
        for(int i=1;i<N-1;i++){
            result(i)=-y(i-1)+2*y(i)-y(i+1)+h*h*(pow(y(i),3)-f(X(i)));
        }
        return result;
    }
    Matrix<double> Jacobian(Vector<double> y){
        Matrix<double> result(N,N);
        //Zeile 1
        result(0,0)=2+h*h*3*y(0)*y(0);
        result(0,1)=-1;
        //Zeile N
        result(N-1,N-2)=-1;
        result(N-1,N-1)=2+h*h*3*y(N-1)*y(N-1);
        //Zeilen 2..N-1
        for(int i=1;i<N-1;i++){
            result(i,i-1)=-1,
            result(i,i)=2+h*h*3*y(i)*y(i);
            result(i,i+1)=-1;
        }
        return result;

    }

    Vector<double> Test(std::function<double(double)> solution){
        Vector<double> result(N);
        for(int i=0;i<N;i++){
            result(i)=solution(X(i));
        }
        return result;
    }

    double X(int i){
        if(i>N-1 || i < 0)
            throw std::invalid_argument("index of X must be in range 0..N-1");
        return h*(i+1);
    }
};

template <typename S>
void Newton(S system, VectorView<double> u, double tol=0.00000001, int max_steps=100){
    Vector<double> u_new = u - inverse(system.Jacobian(u))*system.Eval(u);
    double delta_old = Norm(u_new-u);
    u = u_new;
    for(int i=0; i<max_steps;i++){
        Vector<double> u_new = u - inverse(system.Jacobian(u))*system.Eval(u);
        double delta = Norm(u_new-u);
        u = u_new;
        double q = delta/delta_old;
        double c = q/(1-q)*delta;
        if(q>1) throw std::invalid_argument("newton does not converge");
        if(c<tol) break;
        delta_old = delta;
    }
}


int main(){
    
    std::function sol = [](double x){return sin(M_PI*x);};
    std::function f = [](double x){return M_PI*M_PI*sin(M_PI*x)+pow(sin(M_PI*x),3);};
    FiniteDifference sys(f,50);
    Vector<double> u(50);
    //u=sys.Test(sol);
    Newton(sys,u);
    std::cout << u << std::endl;
    std::cout << sys.Test(sol) << std::endl;
}