#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>
#include <math.h>
#include "expression.h"


namespace ASC_bla {
    template<typename T=double, typename TDIST = std::integral_constant<size_t, 1> >
    class VectorView;
    template<typename T=double>
    class Vector;


    template<typename T, typename TDIST>
    class VectorView : public VecExpr<VectorView<T, TDIST>> {
    protected:
        T *data_;
        size_t size_;
        TDIST dist_;
    public:
        VectorView(size_t size, T *data)
                : data_(data), size_(size) {}

        VectorView(size_t size, TDIST dist, T *data)
                : data_(data), size_(size), dist_(dist) {}

        VectorView(const VectorView& v) // Copy ctor
                : data_(v.data_), size_(v.size_), dist_(v.dist_) {}

        VectorView(const VectorView&& v) // Move ctor
                : data_(v.data_), size_(v.size_), dist_(v.dist_) {}

        VectorView &operator=(const VectorView &v2) { // copy assign
            if(size_ != v2.size_)
                throw std::invalid_argument("Vector size must mazch");
            for (size_t i = 0; i < size_; i++)
                data_[dist_ * i] = v2(i);
            return *this;
        }


        VectorView &operator=(VectorView &&v2) { // move assign
            if(size_ != v2.size_)
                throw std::invalid_argument("Vector size must mazch");
            for (size_t i = 0; i < size_; i++)
                data_[dist_ * i] = v2(i);
            return *this;
        }


        template<typename TB>
        VectorView &operator=(const VecExpr<TB> &v2) {
            for (size_t i = 0; i < size_; i++)
                data_[dist_ * i] = v2(i);
            return *this;
        }

        template<typename TB>
        VectorView &operator+=(const VecExpr<TB> &v2) {
            Vector<T> res{Size()};
            res = (*this)+v2;
            return (*this)=res;
        }
        VectorView &operator*=(T s) {
            Vector<T> res{Size()};
            res = s*(*this);
            return (*this)=res;
        }
        template<typename TB>
        VectorView &operator-=(const VecExpr<TB> &v2) {
            Vector<T> res{Size()};
            res = (*this)+ (-1)*v2;
            return (*this)=res;
        }

        VectorView & operator= (T scal)
        {
          for (size_t i = 0; i < size_; i++)
            data_[dist_*i] = scal;
          return *this;
        }
        T* Data(){return data_;}
        auto View() const { return VectorView(size_, dist_, data_); }
        size_t Size() const { return size_; }
        auto Dist() const { return dist_; }
        T & operator()(size_t i) { return data_[dist_*i]; }
        const T & operator()(size_t i) const { return data_[dist_*i]; }
        
        auto Range(size_t first, size_t next) const {
          return VectorView(next-first, dist_, data_+first*dist_);
        }

        auto Slice(size_t first, size_t slice) const {
            return VectorView<T, size_t>(size_ / slice, dist_ * slice, data_ + first * dist_);
        }

    };


    template<typename T>
    class Vector : public VectorView<T> {
        typedef VectorView<T> BASE;
        using BASE::size_;
        using BASE::data_;
    public:
        explicit Vector(size_t size) 
                : VectorView<T>(size, new T[size]()) { ; }

        Vector(const Vector &v)
                : Vector(v.Size()) {
            *this = v;
        }

        Vector(Vector &&v)
                : VectorView<T>(0, nullptr) {
            std::swap(size_, v.size_);
            std::swap(data_, v.data_);
        }

        template<typename TB>
        Vector(const VecExpr<TB> &v)
                : Vector(v.Size()) {
            *this = v;
        }

        Vector & operator= (T scal)
        {
          for (size_t i = 0; i < size_; i++)
            data_[i] = scal;
          return *this;
        }

        ~Vector() {
             delete[] data_; 
        }

        template<typename TB>
        Vector &operator=(const VecExpr<TB> &v2) {
            for (size_t i = 0; i < size_; i++)
                data_[i] = v2(i);
            return *this;
        }

        Vector &operator=(const Vector &v2) {
            for (size_t i = 0; i < size_; i++)
                data_[i] = v2(i);
            return *this;
        }

        Vector &operator=(Vector &&v2) {
            for (size_t i = 0; i < size_; i++)
                data_[i] = v2(i);
            return *this;
        }

    };


    template<typename ...Args>
    std::ostream &operator<<(std::ostream &ost, const VectorView<Args...> &v) {
        if (v.Size() > 0)
            ost << v(0);
        for (size_t i = 1; i < v.Size(); i++)
            ost << ", " << v(i);
        return ost;
    }
    template<typename T>
    auto Norm(const VecExpr<T>& v){
        decltype(v(0)) sum=0;
        for(int i=0;i<v.Size();i++)
            sum+=v(i)*v(i);
        return sqrt(sum);
        
    }

    template<size_t S, typename T = double>
    class Vec : public VecExpr<Vec<S,T>> {
        T data[S];
    public:
        T & operator()(size_t row) {
            return data[row];
        }

        const T & operator()(size_t row) const {
            return data[row];
        }

        size_t Size() const { return S; }
    };
}

#endif
