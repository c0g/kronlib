#include <functional>
#include <cmath>


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
