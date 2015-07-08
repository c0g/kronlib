//
// Created by Thomas Nickson on 11/06/15.
//

#ifndef KRONMAT_STATIC_FUNCTIONS_H
#define KRONMAT_STATIC_FUNCTIONS_H

#include <vector>


template <typename T>
std::vector<matrix<T>> kronmat_dot_kronmat (
       const std::vector<matrix<T>> & M1, const std::vector<matrix<T>> & M2
)
{
   std::vector<matrix<T>> ans;
   assert(M1.size() == M2.size());
   for (int i = 0; i < M1.size(); ++i) {
       ans.push_back(M1[i] * M2[i]);
   }
   return ans;
}

template <typename T>
matrix<T> kronmat_dot_fullvec ( const std::vector<matrix<T>> & K, const matrix<T> & V )
{
    assert(V.nC() == 1);
    long nmats = K.size();
    matrix<T> x = V;

    for (int n = 0; n < nmats; ++n)
    {
       long thisC = K[n].nC();
       long xSize = x.nR() * x.nC();
       matrix<T> Zt = (K[n] * x.mutable_reshape(thisC, xSize / thisC)).transpose();
       x = Zt.mutable_reshape(Zt.nR() * Zt.nC(), 1);
    }
return x;
}

template <typename T>
matrix<T> kron_full(
       const std::vector<matrix<T>> & K
) {
   long nr = 1;
   long nc = 1;
   long nmats = K.size();
   for (int i = 0; i < nmats; ++i) {
       nr *= K[i].nR();
       nc *= K[i].nC();
   }
   matrix<T> full_mat(nr, nc);
   full_mat = 1;

   long row_acc = 1;
   long col_acc = 1;
   std::deque<long> rowstrides;
   std::deque<long> colstrides;

   rowstrides.push_front(1);
   colstrides.push_front(1);
   for (long n = nmats - 1; n >= 0; --n) {
       row_acc *= K[n].nR();
       col_acc *= K[n].nC();
       rowstrides.push_front(row_acc);
       colstrides.push_front(col_acc);
   }

   for (long rout = 0; rout < nr; ++rout) {
       for (long cout = 0; cout < nc; ++cout) {
           for (long d = 0; d < nmats; ++d) {
               long rowm = (rout % rowstrides[d]) / rowstrides[d+1];
               long colm = (cout % colstrides[d]) / colstrides[d+1];
               full_mat(rout, cout) *= K[d](rowm, colm);
           }
       }
   }
   return full_mat;
}

template <typename T>
matrix<T> kvs_full( const std::vector<matrix<T>> & KVS,
       long start, long end) {
   long nr = end - start;
   long nc = 1;
   long nmats = KVS.size();
   for (int i = 0; i < nmats; ++i) {
       nc *= KVS[i].nC();
   }
   matrix<T> full_mat(nr, nc);
   full_mat = 1;

   long col_acc = 1;
   std::deque<long> colstrides;

   colstrides.push_front(1);
   for (long n = nmats - 1; n >= 0; --n) {
       col_acc *= KVS[n].nC();
       colstrides.push_front(col_acc);
   }

   for (long rout = start; rout < end; ++rout) {
       for (long cout = 0; cout < nc; ++cout) {
           for (long d = 0; d < nmats; ++d) {
               long colm = (cout % colstrides[d]) / colstrides[d+1];
               full_mat(rout - start, cout) *= KVS[d](rout, colm);
           }
       }
   }
   return full_mat;
}

template <typename T>
matrix<T> kvs_full(const std::vector<matrix<T>> & KVS) {
   return kvs_full(KVS, 0, KVS[0].nR());
}
/*
template <typename T>
std::vector<T> vec_dot_kvs(const basic_kvs<T> & kvs, const std::vector<T> & obs)
{
    std::vector<long> shapes = kvs.shapes;
    long full_width = 1;
    for (auto width : shapes) full_width *= width;


    std::vector<T> kroned_vecs;
    kroned_vecs.resize(full_width, 0);


    for (int n = 0; n < kvs.height; ++n) {
        T alpha = -2 * obs[n];
        long current_width = shapes[0];
        std::vector<T> kronme(
            kvs.data[0].begin() + n * shapes[0],
            kvs.data[0].begin() + (n+1) * shapes[0]
        );
        for (int d = 1; d < kvs.dim - 1; ++d) {
            std::vector<T> data(kvs.data[d].begin() + n * shapes[d],kvs.data[d].begin() + (n  + 1) * shapes[d] );
            std::vector<T> destination;
            destination.resize(current_width * shapes[d], 0);

            blas_ger(current_width, shapes[d], alpha,
                     kronme.data(), 1,
                     data.data(), 1,
                     destination.data(), shapes[d]);
            alpha = 1.0;
            kronme = destination;
            current_width *= shapes[d];
        }
        std::vector<T> data (
            kvs.data[kvs.dim-1].begin() + n * shapes[kvs.dim-1],
            kvs.data[kvs.dim-1].begin() + (n  + 1) * shapes[kvs.dim-1]
        );
        blas_ger(current_width, shapes[kvs.dim - 1], 1.0,
                kronme.data(), 1,
                data.data(), 1,
                kroned_vecs.data(), shapes[kvs.dim-1]);
    }
    return kroned_vecs;
}
*/


#endif //KRONMAT_STATIC_FUNCTIONS_H