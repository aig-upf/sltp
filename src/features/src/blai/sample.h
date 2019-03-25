
#pragma once

#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include <common/utils.h>
#include <blai/matrix.h>
#include <blai/transitions.h>
#include "utils.h"

namespace Sample {



class Sample {
public:
    const Matrix matrix_;
    const Transitions transitions_;

    Sample(Matrix&& matrix,
           Transitions&& transitions,
           bool verbose = true)
      : matrix_(matrix), transitions_(transitions)
    {}
    virtual ~Sample() = default;

    const Matrix& matrix() const { return matrix_; }

    const Transitions& transitions() const { return transitions_; }

    bool expanded(unsigned s) const { // CHECK: HACK: this should be in matrix_, like goal()
        return transitions_.num_transitions(s) > 0;
    }

    bool marked(unsigned src, unsigned dst) const {
        return transitions_.marked(src, dst);
    }
};

} // namespaces
