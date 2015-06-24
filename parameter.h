//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_PARAMETER_H
#define KRONMAT_PARAMETER_H

#include <assert.h>
#include <math.h>

enum WARP { POSITIVE, UPPER_BOUND, LOWER_BOUND,  };

template <typename P, typename T>
 P make_param(T value);

template <typename T>
struct warped {
    virtual T val() const = 0;
    virtual T rawval() const = 0;
    virtual void set_val(const T val) = 0;
    virtual void set_rawval(const T val) = 0;
    virtual T dval_drawval() const = 0;
};

template <typename T>
struct positive : warped<T> {
    T raw_value;
    positive() : raw_value(0) {};
    positive(T val) : raw_value(val) {};
    T val() const {
        return exp(raw_value);
    }
    T rawval() const {
        return raw_value;
    }
    void set_rawval(const T value) {
        raw_value = value;
    }
    void set_val(const T value) {
        assert(value > 0);
        raw_value = log(value);
    }
    T dval_drawval() const {
        return exp(raw_value);
    }

};
#endif //KRONMAT_PARAMETER_H
