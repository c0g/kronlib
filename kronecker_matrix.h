////
//// Created by Thomas Nickson on 11/06/15.
////
//
//#ifndef KRONMAT_KRONMAT_H
//#define KRONMAT_KRONMAT_H
//
//#include <vector>
//#include <deque>
//#include <dlib/matrix.h>
//#include "static_functions.h"
//
//
//template <typename T>
//class kronecker_matrix  {
//
//public:
//
//
//private:
//public:
//    kronecker_matrix(const std::vector<dlib::matrix<T>> &sub_matrices) : sub_matrices(sub_matrices) { }
//    kronecker_matrix() {}
//
//public:
//    std::vector<dlib::matrix<T>> sub_matrices;
//    std::vector<long> nRs;
//    std::vector<long> nCs;
//
//
//private:
//    struct literal_assign_helper
//    {
//
//        explicit literal_assign_helper(kronecker_matrix * km_): km(km_) { }
//        ~literal_assign_helper()
//        {
//        }
//
//        const literal_assign_helper& operator, (
//                const dlib::matrix<T> & mat
//        ) const
//        {
//            (*km).push_matrix(mat);
//            return *this;
//        }
//
//    private:
//
//        friend class kronecker_matrix;
//        kronecker_matrix * km;
//
//    };
//
//public:
//
//    kronecker_matrix& operator = (
//            const literal_assign_helper& val
//    )
//    {
//        *this = *val.km;
//        return *this;
//    }
//
//    const literal_assign_helper operator = (
//            const dlib::matrix<T> & mat
//    )
//    {
//        push_matrix(mat);
//        return literal_assign_helper(this);
//    }
//    void push_matrix(const dlib::matrix<T> & mat)
//    {
//        sub_matrices.push_back(mat);
//    }
//
//    kronecker_matrix<T> operator*(
//            const kronecker_matrix<T> & other
//    )
//    {
//        kronecker_matrix ans;
//        ans.sub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
//        return ans;
//    }
//    long nr()
//    {
//        long rows = 1;
//        for (auto m : sub_matrices)
//        {
//            rows *= m.nr();
//        }
//        return rows;
//    }
//    long nc()
//    {
//        long cols = 1;
//        for (auto m : sub_matrices)
//        {
//            cols *= m.nc();
//        }
//        return cols;
//    }
//    dlib::matrix<T> full() {
//        dlib::matrix<T> full_mat = dlib::ones_matrix<T>(nr(), nc());
//        long row_acc = 1;
//        long col_acc = 1;
//        long nmats = sub_matrices.size();
//        std::deque<long> rowstrides;
//        std::deque<long> colstrides;
//
//        rowstrides.push_front(1);
//        colstrides.push_front(1);
//        for (long n = nmats - 1; n >= 0; --n) {
//            row_acc *= sub_matrices[n].nr();
//            col_acc *= sub_matrices[n].nc();
//            rowstrides.push_front(row_acc);
//            colstrides.push_front(col_acc);
//        }
//
//        for (long rout = 0; rout < nr(); ++rout) {
//            for (long cout = 0; cout < nc(); ++cout) {
//                for (long d = 0; d < nmats; ++d) {
//                    long rowm = (rout % rowstrides[d]) / rowstrides[d+1];
//                    long colm = (cout % colstrides[d]) / colstrides[d+1];
//                    full_mat(rout, cout) *= sub_matrices[d](rowm, colm);
//                }
//            }
//        }
//        return full_mat;
//    }
//
//};
//template <typename T>
//std::ostream& operator<<(
//        std::ostream& out, kronecker_matrix<T> K
//)
//{
//    int matnum = 1;
//    for (auto m : K.sub_matrices) {
//        out << "Submatrix: " << matnum++ << std::endl;
//        out << m << std::endl;
//    }
//    return out;
//}
//
//#endif //KRONMAT_KRONMAT_H
