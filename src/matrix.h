#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>
#include <vector>

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
        Matrix(size_t width, size_t height) : _data(width * height), _width{width}, _height{height} {}

        T &operator()(size_t row, size_t column) {
            if (ORD == Ordering::ColMajor) {
                return _data[row * _height + column];
            } else {
                return _data[column * _width + row];
            }
        }

        int Width() { return _width; }

        int Height() { return _height; }

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

    }


    template<typename T>
    Vector <T> operator*(const Matrix<T> &a, const Vector <T> &x) {

    }

    template<typename T>
    Vector <T> operator*(const Vector <T> &x, const Matrix<T> &a) {

    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &a, T c) {

    }

}
#endif