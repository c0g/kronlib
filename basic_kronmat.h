//
// Created by Thomas Nickson on 15/06/15.
//

#ifndef KRONMAT_BASIC_KRONMAT_H
#define KRONMAT_BASIC_KRONMAT_H

#include <vector>

template <typename T = float>
struct basic_kvs {
public:
    long height;
    long dim;
    std::vector<long> shapes;
    std::vector<std::vector<T>> data;


public:
    basic_kvs<T>(long height_, std::vector<long> shapes_, std::vector<std::vector<T>> data_) :
            height(height_), shapes(shapes_), data(data_), dim(shapes_.size()) {}

};


#endif //KRONMAT_BASIC_KRONMAT_H
