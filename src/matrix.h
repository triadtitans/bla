#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>
#include <vector>
#include <exception>
#include <utility>

#include "expression.h"
#include "vector.h"

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
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    MatrixView &operator=(T scal) {
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= scal;
            }
        }
        return *this;
    }

    auto View() const { return MatrixView(_height,_width,_data, _dist); }

    size_t Height() const { return _height; }
    size_t Width() const { return _width; }
    size_t Dist() const { return _dist; }

    MatrixView Cols(size_t first, size_t next) {
        // constructor: size_t height, size_t width, T* data, size_t dist
        return MatrixView(
            this->Height(),
            next,
            this->_data + first,
            this->_dist
        );
    }

    MatrixView Rows(size_t first, size_t next) {
        return MatrixView(
            next,
            this->Width(),
            this->_data + this->_dist * first,
            this->_dist
        );
    }

    VectorView<T, size_t> Row(size_t i) {
        if(ORD == Ordering::RowMajor) {
            return VectorView<T, size_t>(_width, 1, _data+i*_width);
        }

        return VectorView<T, size_t>(_width, _height, _data+i);
    }

    VectorView<T, size_t> Col(size_t j) {
        if(ORD == Ordering::ColMajor) {
            return VectorView<T, size_t>(_height, 1, _data+j*_width);
        }

        return VectorView<T, size_t>(_height, _width, _data+j);
    }

    MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor> Transpose() {
        return MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor>(_height, _width, _data, ORD == Ordering::RowMajor ? _width : _height);
    };
    void RowSwap(size_t row,size_t dest){
        for(size_t i=0; i<Width();i++)
            std::swap((*this)(dest,i),(*this)(row,i));
    }
    void EnsureNonzero(size_t row, size_t column){
        size_t i=row;
        for(; i<Height();i++)
            if((*this)(i,column)!=0) break;
        if(i==Height())
            throw std::invalid_argument("Matrix singular");
        RowSwap(column,i);
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
    using BASE::_dist;

public:
    Matrix(size_t height, size_t width)
            : MatrixView<T,ORD>(height, width, new T[height*width], ORD == Ordering::RowMajor ? width : height ) {}

    Matrix(const Matrix &m)
            : Matrix(m.Height(),m.Width()) {
        *this = m;
    }

    Matrix(Matrix &&m)
            : MatrixView<T,ORD>{0,0, nullptr,0} {
        std::swap(_width, m._width);
        std::swap(_height, m._height);
        std::swap(_data, m._data);
        std::swap(_dist,m._dist);
    }

    template<typename TA>
    Matrix(const MatrixExpr<TA> &m)
            : Matrix(m.Height(),m.Width()) {
        *this = m;
    }


    ~Matrix() { delete[] _data; }

    using BASE::operator=;

    Matrix &operator=(const Matrix &m2) {
        _width = this->_width;
        _height = this->_height;
        _dist = this->_dist;
        for (size_t i = 0; i < _width*_height; i++){
            _data[i]=m2._data[i];
        }
        return *this;
    }

    Matrix &operator=(Matrix &&m) {
        _width = 0;
        _height = 0;
        _data = nullptr;
        _dist = 0;
        std::swap(_height, m._height);
        std::swap(_width, m._width);
        std::swap(_data, m._data);
        std::swap(_dist, m._dist);
        return *this;
    }
};

template<typename T, Ordering ORD>
MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor> Transpose(MatrixView<T, ORD> &m) {
    return m.Transpose();
}


template <typename T>
Matrix<T> inverse(const MatrixView<T>& m ){
    if(m.Height() != m.Width())
        throw std::invalid_argument("Only square matrixes allowed");
    Matrix<T> work{m.Height(),m.Height()*2};
    for(size_t row=0; row < work.Height();row++){
        for(size_t column=0; column < work.Height();column++){
            work(row,column)=m(row,column);
            work(row,column+m.Height()) = (row==column) ? 1 : 0;
        }
    }
    std::cout << work;
    for(size_t column=0; column < work.Height();column++){
        work.EnsureNonzero(column,column);
        std::cout << work;
        work.Rows(column,column+1)=1.0/(work(column,column))*work.Rows(column,column+1);
        //work.RowMultiply(column,1/(work(column,column)));
        std::cout << work;
        for(size_t row=0; row < work.Height();row++){
            if(row==column) continue;
            work.Rows(row, row+1) = work.Rows(row, row+1) + (-work(row,column))*work.Rows(column,column+1);
            //work = work.RowMulAdd(column,row,-work(row,column));
            std::cout << work;
        }
    }
    std::cout << work;
    Matrix<T> result = work.Cols(m.Width(),m.Width());
    return result;
}

}
#endif