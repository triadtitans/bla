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

        T &operator()(size_t row, size_t column) {
            if constexpr (ORD == Ordering::ColMajor) {
                return _data[row * _height + column];
            } else {
                return _data[column * _width + row];
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

    }

    template<typename T>
    Vector<T> operator*(const Vector<T> &x, const Matrix<T> &a) {

    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &a, T c) {

    }

}
#endif