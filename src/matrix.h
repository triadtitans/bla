#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>
#include <vector>
#include <exception>

namespace ASC_bla {
    enum class Ordering {
        ColMajor, RowMajor
    };

    template<typename T, Ordering ORD = Ordering::RowMajor>
    class Matrix {
        std::vector<T> _data;
        size_t _width;
        size_t _height;
    public:
        Matrix(size_t height,size_t width) : _data(width * height), _width{width}, _height{height} {}
        Matrix(const Vector<T>& v) : _data(v.Size()), _width{1}, _height{v.Size()} {
            for(size_t i=0;i<v.Size();i++)
                (*this)(i,0)=v(i);
        }

        T& operator()(size_t row, size_t column) {
            if constexpr (ORD == Ordering::ColMajor) {
                return _data[column * _height + row];
            } else {
                return _data[row * _width + column];
            }
        }

        const T& operator()(size_t row, size_t column) const {
            if constexpr (ORD == Ordering::ColMajor) {
                return _data[row * _height + column];
            } else {
                return _data[column * _width + row];
            }
        }

        size_t Width() const { return _width; }

        size_t Height() const { return _height; }

        Matrix<T> Inverse() const{
            if(this->Height() != this->Width())
                throw std::invalid_argument("Only square matrixes allowed");
            Matrix<T> work{this->Height(),this->Height()*2};
            for(size_t row=0; row < work.Height();row++){
                for(size_t column=0; column < work.Height();column++){
                    work(row,column)=(*this)(row,column);
                    work(row,column+Height()) = (row==column) ? 1 : 0;
                }
            }
            work.Dump();
            return work;
        }

        void Dump() {
            for (size_t r = 0; r < this->Height(); r++) {
                for (size_t c = 0; c < this->Width(); c++) {
                    std::cout << (*this)(r, c) << " ";
                }
                std::cout << " \n";
            }
        }
    };

    template<typename T>
    Matrix<T> operator+(const Matrix<T> &a, const Matrix<T> &b) {

    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b) {
        if(a.Width() != b.Height()){
            throw std::invalid_argument("Matrix dimension dont match");
        }
        Matrix<T> result{a.Height(),b.Width()};
        for(size_t row=0; row < result.Height();row++){
            for(size_t column=0; column < result.Width();column++){
                //Calculate C_row,column
                T sum=0; //TODO: Zero element of T
                for(size_t i=0;i<a.Width();i++){
                    sum += a(row,i)*b(i,column);
                }
                result(row,column)=sum;
            } 
        }
        return result;
    }


    template<typename T>
    Vector<T> operator*(const Matrix<T> &a, const Vector<T> &x) {
        Matrix m = a*Matrix(x);
        Vector<T> res{m.Height()};
        for(size_t i=0;i<res.Size();i++)
            res(i)=m(i,0);
        return res;
    }

    template<typename T>
    Vector<T> operator*(const Vector<T> &x, const Matrix<T> &a) {

    }

    template<typename T>
    Matrix<T> operator*(T c, const Matrix<T> &a) {
        Matrix<T> result = a;
        for(size_t row=0; row < result.Height();row++){
            for(size_t column=0; column < result.Width();column++){
                result(row,column)*=c;
            } 
        }
        return result;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &a, T c) {
        return c*a;
    }

}
#endif