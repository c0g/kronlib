////
//// Created by Thomas Nickson on 11/06/15.
////
//
//#ifndef KRONMAT_STATIC_FUNCTIONS_H
//#define KRONMAT_STATIC_FUNCTIONS_H
//
//#include <vector>
//#include <dlib/matrix.h>
//
//
//template <typename T>
//std::vector<dlib::matrix<T>> kronmat_dot_kronmat (
//        const std::vector<dlib::matrix<T>> & M1, const std::vector<dlib::matrix<T>> & M2
//)
//{
//    std::vector<dlib::matrix<T>> ans;
//    DLIB_CASSERT(M1.size() == M2.size(), "M1 and M2 must be the same size");
//    for (int i = 0; i < M1.size(); ++i) {
//        ans.push_back(M1[i] * M2[i]);
//    }
//    return ans;
//}
//
//template <typename T>
//dlib::matrix<T> kronmat_dot_fullvec (
//        const std::vector<dlib::matrix<T>> & K, const dlib::matrix<T> & V
//)
//{
//    long nmats = K.size();
//    dlib::matrix<T> x = V;
//    for (int n = 0; n < nmats; ++n)
//    {
//        long thisC = K[n].nc();
//        long xSize = x.nr() * x.nc();
//        dlib::matrix<T> Zt = dlib::trans(K[n] * dlib::reshape(x, thisC, xSize / thisC));
//        x = dlib::reshape_to_column_vector(Zt);
//    }
//    return x;
//}
//
//template <typename T>
//dlib::matrix<T> kron_full(
//        const std::vector<dlib::matrix<T>> & K
//) {
//    long nr = 1;
//    long nc = 1;
//    long nmats = K.size();
//    for (int i = 0; i < nmats; ++i) {
//        nr *= K[i].nr();
//        nc *= K[i].nc();
//    }
//    dlib::matrix<T> full_mat = dlib::ones_matrix<T>(nr, nc);
//    long row_acc = 1;
//    long col_acc = 1;
//    std::deque<long> rowstrides;
//    std::deque<long> colstrides;
//
//    rowstrides.push_front(1);
//    colstrides.push_front(1);
//    for (long n = nmats - 1; n >= 0; --n) {
//        row_acc *= K[n].nr();
//        col_acc *= K[n].nc();
//        rowstrides.push_front(row_acc);
//        colstrides.push_front(col_acc);
//    }
//
//    for (long rout = 0; rout < nr; ++rout) {
//        for (long cout = 0; cout < nc; ++cout) {
//            for (long d = 0; d < nmats; ++d) {
//                long rowm = (rout % rowstrides[d]) / rowstrides[d+1];
//                long colm = (cout % colstrides[d]) / colstrides[d+1];
//                full_mat(rout, cout) *= K[d](rowm, colm);
//            }
//        }
//    }
//    return full_mat;
//}
//
//template <typename T>
//dlib::matrix<T> kvs_full(
//        const std::vector<dlib::matrix<T>> & KVS,
//        long start, long end
//) {
//    long nr = end - start;
//    long nc = 1;
//    long nmats = KVS.size();
//    for (int i = 0; i < nmats; ++i) {
//        nc *= KVS[i].nc();
//    }
//    dlib::matrix<T> full_mat = dlib::ones_matrix<T>(nr, nc);
//
//    long col_acc = 1;
//    std::deque<long> colstrides;
//
//    colstrides.push_front(1);
//    for (long n = nmats - 1; n >= 0; --n) {
//        col_acc *= KVS[n].nc();
//        colstrides.push_front(col_acc);
//    }
//
//    for (long rout = start; rout < end; ++rout) {
//        for (long cout = 0; cout < nc; ++cout) {
//            for (long d = 0; d < nmats; ++d) {
//                long colm = (cout % colstrides[d]) / colstrides[d+1];
//                full_mat(rout - start, cout) *= KVS[d](rout, colm);
//            }
//        }
//    }
//    return full_mat;
//}
//
//template <typename T>
//dlib::matrix<T> kvs_full(
//        const std::vector<dlib::matrix<T>> & KVS
//) {
//    return kvs_full(KVS, 0, KVS[0].nr());
//}
//
//template <typename T>
//dlib::matrix<T> kronvecstack_dot_fullvec (
//        const std::vector<dlib::matrix<T>> & KVS, const dlib::matrix<T> & V
//)
//{
//
//}
//
//#endif //KRONMAT_STATIC_FUNCTIONS_H