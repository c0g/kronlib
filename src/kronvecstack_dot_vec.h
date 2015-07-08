// //
// // Created by Thomas Nickson on 15/06/15.
// //
// #ifndef KRONMAT_KRONVECSTACK_DOT_VEC_H
// #define KRONMAT_KRONVECSTACK_DOT_VEC_H

// #include <vector>
// #include <iostream>
// #include <deque>
// // #include "basic_kronmat.h"
// #include "blas_wrap.h"


// template <typename T>
// std::vector<T> kronvecstack_dot_vec(const basic_kvs<T> & kvs, const std::vector<T> & vec) {
//     // Size of initial reshape of vec
//     long init_width = 1;
//     long last_idx = kvs.dim - 1;
//     for (int i = 0; i < last_idx; ++i) {
//         init_width *= kvs.shapes[i];
//     }

//     std::vector<T> alpha;
//     alpha.resize(kvs.height * init_width);

// //    We can use BLAS for the first part. Gives us a N x \prod_{i=1}^{D-1} dim_i matrix.
//     blas_gemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//               kvs.height, init_width, kvs.shapes[last_idx],
//               1.0,
//               kvs.data[last_idx].data(), kvs.shapes[last_idx],
//               vec.data(), kvs.shapes[last_idx],
//               1.0,
//               alpha.data(), init_width);

//     std::vector<float> beta;
//     long next_width = init_width / kvs.shapes[last_idx];
//     long next_idx = last_idx - 1;
//     beta.resize(kvs.height * next_width);

//     long N = kvs.height;
//     std::vector<long> shapes = kvs.shapes;
//     long full_width = 1;
//     for (auto width : shapes) full_width *= width;
//     long current_width = full_width / shapes.back();


//     for (long d = kvs.dim - 2; d >= 0; --d)
//     {
//         current_width = current_width / shapes[d];
//         std::vector<T> tmp;
//         int this_width = shapes[d];

//         for (long n = 0; n < kvs.height; ++n)
//         {
//             for (long idx = 0; idx < current_width; ++idx)
//             {
//                 T val = blas_dot(this_width,
//                             &kvs.data[d][n * this_width], 1,
//                             &alpha[n * current_width * this_width + idx * this_width ], 1);
//                 tmp.push_back(val);
//             }
//         }
//         alpha = tmp;
//     }
//     return alpha;
// }
// #endif //KRONMAT_KRONVECSTACK_DOT_VEC_H
