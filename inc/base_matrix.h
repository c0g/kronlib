template <typename T>
class base_matrix {

    public:
        virtual long nR() const = 0;
        virtual long nC() const = 0;
        virtual bool isTrans() const = 0;
        virtual T trace() const = 0;
        virtual base_matrix<T> transpose() const = 0;
        virtual void mutable_transpose()  = 0;
        virtual T& operator()(long ridx, long cidx) = 0;
        virtual T operator()(long ridx, long cidx) const  = 0;
        virtual base_matrix<T> operator*(const base_matrix<T> &other) const = 0;
        virtual base_matrix<T> Tdot(const base_matrix<T> &other)  = 0;
        virtual base_matrix<T> hadamard(const base_matrix<T> &other) const = 0;
        virtual base_matrix<T> elemwise_div(const base_matrix<T> &other) const = 0;
        virtual void operator*=(const T val) = 0;
        virtual void operator/=(const T val)  = 0;
        virtual base_matrix<T> operator*(const T val) const = 0;
        virtual base_matrix<T> operator/(const T val) const = 0;
        virtual base_matrix<T> operator==(const base_matrix<T> & other) const = 0;
        virtual bool operator!=(const base_matrix<T> & other) const = 0;
        virtual void exp_inplace() = 0;
    };
