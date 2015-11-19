#include <functional>
#include <cmath>
#include <limits>


template <typename T>
struct tol_equal : public std::binary_function<T, T, bool>
{

    T tol;
    __host__ __device__
        tol_equal(T tol) : tol{tol} {}
    __host__ __device__
    bool operator()(T lhs, T rhs)
    {
        return std::abs(lhs - rhs) <= tol;
    }
};

template <typename T>
struct div_by : public std::unary_function<T, T>
{
    T num;
    __host__ __device__
    div_by(T _num) : num{_num} {};

    __host__ __device__
    T operator()(T val)
    {
        return val / num;
    }
};

template <typename T>
struct add_by : public std::unary_function<T, T>
{
    T num;
    __host__ __device__
    add_by(T _num) : num{_num} {};

    __host__ __device__
    T operator()(T val)
    {
        return val + num;
    }
};

template <typename T>
struct minus_by : public std::unary_function<T, T>
{
    T num;
    __host__ __device__
    minus_by(T _num) : num{_num} {};

    __host__ __device__
    T operator()(T val)
    {
        return val - num;
    }
};

template <typename T>
struct mult_by : public std::unary_function<T, T>
{
    T num;
    __host__ __device__
    mult_by(T _num) : num{_num} {};

    __host__ __device__
    T operator()(T val)
    {
        return val * num;
    }
};

template <typename T>
struct logarithm : public std::unary_function<T, T>
{

    __host__ __device__
    T operator()(T val)
    {
        return log(val);
    }
};

template <typename T>
struct exponentiate : public std::unary_function<T, T>
{

    __host__ __device__
    T operator()(T val)
    {
        return exp(val);
    }
};

struct and_reduce : public std::binary_function<bool,bool,bool>
{

    __host__ __device__
    bool operator()(bool left, bool right)
    {
        return left && right;
    }
};
struct transpose_index : public std::unary_function<size_t,size_t>
{
    size_t nc, nr;

    __host__ __device__
    transpose_index(size_t _nc, size_t _nr) : nc{_nc}, nr{_nr} {}

    __host__ __device__
    size_t operator()(size_t linear_index)
    {
        size_t i = linear_index / nr;
        size_t j = linear_index % nr;

        return nc * j + i;
    }
};

struct kronecker_index : public std::unary_function<size_t, size_t>
{
    size_t out_rs, out_cs, source_rs, source_cs, this_rowstride, next_rowstride, this_colstride, next_colstride;
    __host__ __device__
    kronecker_index(size_t out_rs, size_t out_cs, size_t source_rs, size_t source_cs, 
            size_t this_rowstride, size_t next_rowstride, size_t this_colstride, size_t next_colstride) :
                    out_rs{out_rs}, out_cs{out_cs}, source_rs{source_rs}, source_cs{source_cs}, 
                    this_rowstride{this_rowstride}, next_rowstride{next_rowstride}, this_colstride{this_colstride}, next_colstride{next_colstride} {}

    __host__ __device__ 
    size_t operator()(size_t output_index) 
    {
       size_t rout = output_index % out_rs;
       size_t cout = output_index / out_rs;
       size_t rsource = (rout % this_rowstride) / next_rowstride;
       size_t csource = (cout % this_colstride) / next_colstride;
       return csource * source_rs + rsource;
    }
};
template <typename T>
struct zero_if_not_lower_triangular : public std::binary_function<size_t, T, T>
{
    size_t n;
    __host__ __device__
    zero_if_not_lower_triangular(size_t n) : n{n} {}

    __host__ __device__ 
    T operator()(size_t idx, T val) 
    {
       size_t r = idx % n;;
       size_t c = idx / n;
       T ans = (c <= r) ? val : 0;
       //std::cout << r << ", " << c << ", " << ans << std::endl;
       return ans;
    }
};
