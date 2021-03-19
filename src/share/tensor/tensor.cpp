#ifndef __TENSOR_CPP
#define __TENSOR_CPP

#include "tensor.h"

namespace HyperFlow {

/* Constructor */
Tensor1D::Tensor1D()
{}

/* Constructor with supplied extent */
Tensor1D::Tensor1D(const unsigned int size)
{
    vec = std::vector<double>(size);
}

/* Constructor with supplied extent and default values */
Tensor1D::Tensor1D(const unsigned int size, const double default_val)
{
    vec = std::vector<double>(size, default_val);
}

/* Constructor via initialiser list */
Tensor1D::Tensor1D(std::initializer_list<double> init_vals) : vec(init_vals)
{}

/* Destructor */
Tensor1D::~Tensor1D()
{}

/* Underlying write access to the std::vector container */
std::vector<double>& Tensor1D::get_vec()
{
    return vec;
}

/* Underlying read access to the std::vector container */
const std::vector<double>& Tensor1D::get_vec() const
{
    return vec;
}

/* Element-wise write access to the std::vector container */
double& Tensor1D::operator[](const unsigned int index)
{
    return vec[index];
}

/* Element-wise read access to the std::vector container */
const double& Tensor1D::operator[](const unsigned int index) const
{
    return vec[index];
}

/* Add a new element to the end of the tensor */
void Tensor1D::push_back(const double s)
{
    vec.push_back(s);
}

/* Read access to the size of the tensor */
const unsigned int Tensor1D::size() const
{
    return vec.size();
}
        
/* Inner product with another tensor */
double Tensor1D::dot(const Tensor1D& t) const
{
    double inner_product = 0.0;
    for (unsigned int i=0; i<vec.size(); i++) {
        inner_product += vec[i] * t[i];
    }
    return inner_product;
};

/* L2 (Euclidean) Norm for tensor 'length' */
double Tensor1D::l2_norm() const
{
    return sqrt(this->dot(*this));
}

/* Return a normalised unit tensor in the same direction */
Tensor1D Tensor1D::unit()
{
    double inv_length = 1.0 / this->l2_norm();
    Tensor1D tcpy = (*this) * inv_length;
    return tcpy;
}

/* Console readable output of the tensor */
void Tensor1D::output() const
{
    std::cout << std::setprecision(8) << "[";

    for (unsigned int i=0; i<vec.size(); i++) {
        if (fabs(vec[i]) < 1e-12) {
            std::cout << 0.0;
        } else {
            std::cout << vec[i];
        }
        if (i < vec.size() - 1) {
            std::cout << ", ";    
        }
    }
    
    std::cout << "]" << std::endl;
}

/* Element-wise vector addition */
Tensor1D operator+(const Tensor1D& a, const Tensor1D& b)
{
    Tensor1D result;
    result.get_vec().reserve(a.size());

    std::transform(
        a.get_vec().begin(),
        a.get_vec().end(),
        b.get_vec().begin(),
        std::back_inserter(result.get_vec()), std::plus<double>()
    );

    return result;    
};

/* Element-wise vector incremental addition */
Tensor1D& operator+=(Tensor1D& a, const Tensor1D& b)
{
    for (unsigned int i=0; i<a.size(); i++) {
        a[i] += b[i];
    }
    return a;
};

/* Element-wise vector multiplication */
Tensor1D operator*(const Tensor1D& a, const Tensor1D& b)
{
    Tensor1D result;
    result.get_vec().reserve(a.size());

    std::transform(
        a.get_vec().begin(),
        a.get_vec().end(),
        b.get_vec().begin(),
        std::back_inserter(result.get_vec()), std::multiplies<double>()
    );

    return result;    
};

/* Scalar-vector addition (scalar LHS) */
Tensor1D operator+(const double s, const Tensor1D& t)
{
    Tensor1D tcpy(t.size());

    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        tcpy.get_vec().begin(),
        std::bind1st(std::plus<double>(), s)
    );

    return tcpy;
};

/* Scalar-vector addition (scalar RHS) */
Tensor1D operator+(const Tensor1D& t, const double s)
{
    return s + t;
};

/* Scalar-vector incremental addition (scalar RHS) */
Tensor1D& operator+=(Tensor1D& t, const double s)
{
    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        t.get_vec().begin(),
        std::bind1st(std::plus<double>(), s)
    );

    return t;
};

/* Element-wise vector addition (a - b) */
Tensor1D operator-(const Tensor1D& a, const Tensor1D& b)
{
    Tensor1D result;
    result.get_vec().reserve(a.size());

    std::transform(
        a.get_vec().begin(),
        a.get_vec().end(),
        b.get_vec().begin(),
        std::back_inserter(result.get_vec()), std::minus<double>()
    );

    return result;
};

/* Element-wise vector incremental addition (a -= b) */
Tensor1D& operator-=(Tensor1D& a, const Tensor1D& b)
{
    for (unsigned int i=0; i<a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
};

/* Scalar-vector subtraction (scalar LHS) */
Tensor1D operator-(const double s, const Tensor1D& t)
{
    Tensor1D tcpy(t.size());

    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        tcpy.get_vec().begin(),
        std::bind1st(std::minus<double>(), s)
    );

    return tcpy;
};

/* Scalar-vector subtraction (scalar RHS) */
Tensor1D operator-(const Tensor1D& t, const double s)
{
    return -1.0 * s + t;
};

/* Scalar-vector incremental subtraction (scalar RHS) */
Tensor1D& operator-=(Tensor1D& t, const double s)
{
    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        t.get_vec().begin(),
        std::bind2nd(std::minus<double>(), s)
    );

    return t;
};

/* Scalar-vector multiplication (scalar LHS) */
Tensor1D operator*(const double s, const Tensor1D& t)
{
    Tensor1D tcpy(t.size());

    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        tcpy.get_vec().begin(),
        std::bind1st(std::multiplies<double>(), s)
    );

    return tcpy;
};

/* Scalar-vector multiplication (scalar RHS) */
Tensor1D operator*(const Tensor1D& t, const double s)
{
    return s * t;
};

/* Scalar-vector incremental multiplication (scalar RHS) */
Tensor1D& operator*=(Tensor1D& t, const double s)
{
    std::transform(
        t.get_vec().begin(),
        t.get_vec().end(),
        t.get_vec().begin(),
        std::bind1st(std::multiplies<double>(), s)
    );

    return t;
};

}

#endif
