#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>
#include <vector>
#include <exception>

#include "expression.h"

namespace ASC_bla {
    enum class Ordering {
        ColMajor, RowMajor
    };

template <typename T, Ordering ORD = Ordering::RowMajor>
class MatrixView : public MatrixExpr<MatrixView<T,ORD>>{
protected:
    T* _data;
    size_t _width;
    size_t _height;
    size_t _dist;
public:
    MatrixView(size_t height, size_t width, T* data, size_t dist):
        _data{data},_width{width},_height{height},_dist{dist} {};

    template<typename TB>
    MatrixView &operator=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() && _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                this(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    MatrixView &operator=(T scal) {
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                this(row, col)= scal;
            }
        }
        return *this;
    }

    auto View() const { return MatrixView(_height,_width,_data, _dist); }

    size_t Height() const { return _height; }
    size_t Width() const { return _width; }

    MatrixView Cols(size_t first, size_t next) {
        return MatrixView(
            next,
            this->Height(),
            this->_data + first,
            this->_dist
        );
    }

    MatrixView Rows(size_t first, size_t next) {
        return MatrixView(
            this->Width(),
            next,
            this->_data + this->_dist * first,
            this->_dist
        );
    }

    T& operator()(size_t row, size_t column) {
        return const_cast<T&>(std::as_const(*this)(row,column));
    }

    const T& operator()(size_t row, size_t column) const {
        if constexpr (ORD == Ordering::RowMajor) {
            return _data[row * _dist + column];
        } else {
            return _data[column * _dist + row];
        }
    }
};

template<typename T, Ordering ORD=Ordering::RowMajor>
class Matrix : public MatrixView<T,ORD> {
    typedef MatrixView<T> BASE;
    using BASE::_width;
    using BASE::_height;
    using BASE::_data;

public:
    Matrix(size_t height, size_t width)
            : MatrixView<T,ORD>(height, width, new T[height*width], ORD == Ordering::RowMajor ? width : height ) {}

    Matrix(const Matrix &m)
            : Matrix(m.Height(),m.Width()) {
        *this = m;
    }

    Matrix(Matrix &&m)
            : MatrixView<T,ORD>{0, nullptr} {
        std::swap(_width, m._width);
        std::swap(_height, m._height);
        std::swap(_data, m._data);
    }

    template<typename TA>
    Matrix(const MatrixExpr<TA> &m)
            : Matrix(m.Height()+m.Width()) {
        *this = m;
    }


    ~Matrix() { delete[] _data; }

    using BASE::operator=;

    Matrix &operator=(const Matrix &m2) {
        _width = this->_width;
        _height = this->_height;
        for (size_t i = 0; i < _width*_height; i++){
            _data[i]=m2->_data[i];
        }
        return *this;
    }

    Matrix &operator=(Matrix &&m) {
        _width = 0;
        _height = 0;
        _data = nullptr;
        std::swap(_height, m._height);
        std::swap(_width, m._width);
        std::swap(_data, m._data);
        return *this;
    }
};

}


#endif