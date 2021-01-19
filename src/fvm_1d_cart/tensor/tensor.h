#ifndef __TENSOR_H
#define __TENSOR_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace HyperFlow {

class Tensor1D {

    public:
        
        /* Constructor */
        Tensor1D();

        /* Constructor with supplied extent */
        Tensor1D(const unsigned int size);

        /* Constructor with supplied extent and default values */
        Tensor1D(const unsigned int size, const double default_val);

        /* Constructor via initialiser list */
        Tensor1D(std::initializer_list<double> init_vals);

        /* Destructor */
        virtual ~Tensor1D();

        /* Underlying read/write access to the std::vector container */
        std::vector<double>& get_vec();

        /* Underlying read access to the std::vector container */
        const std::vector<double>& get_vec() const;

        /* Element-wise write access to the std::vector container */
        double& operator[](const unsigned int index);

        /* Element-wise read access to the std::vector container */
        const double& operator[](const unsigned int index) const;

        /* Add a new element to the end of the tensor */
        void push_back(const double s);

        /* Read access to the size of the tensor */
        const unsigned int size() const;
        
        /* Inner product with another tensor */
        double dot(const Tensor1D& t) const;

        /* L2 (Euclidean) Norm for tensor 'length' */
        double l2_norm() const;

        /* Return a normalised unit tensor in the same direction */
        Tensor1D unit();

        /* Console readable output of the tensor */
        void output() const;

    private:
        
        /* Underlying storage vector */
        std::vector<double> vec;
};

/* Element-wise vector addition */
Tensor1D operator+(const Tensor1D& a, const Tensor1D& b);

/* Element-wise vector incremental addition */
Tensor1D& operator+=(Tensor1D& a, const Tensor1D& b);

/* Element-wise vector multiplication */
Tensor1D operator*(const Tensor1D& a, const Tensor1D& b);

/* Scalar-vector addition (scalar LHS) */
Tensor1D operator+(const double s, const Tensor1D& t);

/* Scalar-vector addition (scalar RHS) */
Tensor1D operator+(const Tensor1D& t, const double s);

/* Scalar-vector incremental addition (scalar RHS) */
Tensor1D& operator+=(Tensor1D& t, const double s);

/* Element-wise vector subtraction (a - b) */
Tensor1D operator-(const Tensor1D& a, const Tensor1D& b);

/* Element-wise vector incremental subtraction (a -= b) */
Tensor1D& operator-=(Tensor1D& a, const Tensor1D& b);

/* Scalar-vector subtraction (scalar LHS) */
Tensor1D operator-(const double s, const Tensor1D& t);

/* Scalar-vector subtraction (scalar RHS) */
Tensor1D operator-(const Tensor1D& t, const double s);

/* Scalar-vector incremental subtraction (scalar RHS) */
Tensor1D& operator-=(Tensor1D& t, const double s);

/* Scalar-vector multiplication (scalar LHS) */
Tensor1D operator*(const double s, const Tensor1D& t);

/* Scalar-vector multiplication (scalar RHS) */
Tensor1D operator*(const Tensor1D& t, const double s);

/* Scalar-vector incremental multiplication (scalar RHS) */
Tensor1D& operator*=(Tensor1D& t, const double s);

typedef Tensor1D Vec1D;
typedef std::vector<Tensor1D> Vec2D;
typedef std::vector<std::vector<Tensor1D> > Vec3D;
typedef std::vector<std::vector<std::vector<Tensor1D> > > Vec4D;

}

#endif
