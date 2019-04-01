
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

//! A simple container for a pair of feature matrix and transition sample
class Sample {
public:
    const FeatureMatrix matrix_;
    const TransitionSample transitions_;

    Sample(FeatureMatrix&& matrix, TransitionSample&& transitions)
      : matrix_(matrix), transitions_(transitions)
    {}
    virtual ~Sample() = default;

    const FeatureMatrix& matrix() const { return matrix_; }
    const TransitionSample& transitions() const { return transitions_; }

    //! Remap the states in the current sample, whose IDs are assumed to span the full range [0..n],
    //! into a smaller range [0..m] with the m states that are either contained in `selected` or a successor of them
    Sample* resample(const std::unordered_set<unsigned>& selected) const {

        // Add all successors of selected states
        std::set<unsigned> closed; // set must be ordered
        for (unsigned s:selected) {
            closed.insert(s);
            const auto& succ = transitions_.successors(s);
            closed.insert(succ.begin(), succ.end());
        }

        // Generate the mapping
        std::unordered_map<unsigned, unsigned> mapping;
        unsigned i = 0;
        for (unsigned s:closed) { // closed is guaranteed to be sorted
            mapping.emplace(s, i++);
        }

//        std::cout << "Selected: " << std::endl;
//        for (unsigned s:selected)  std::cout << s << ", ";
//        std::cout << std::endl;
//
//        std::cout << "Closed: " << std::endl;
//        for (unsigned s:closed)  std::cout << s << ", ";
//        std::cout << std::endl;
//
//        std::cout << "Mapping: " << std::endl;
//        for (auto& elem:mapping)  std::cout << elem.first << ": " << elem.second << ", ";
//        std::cout << std::endl;

        return new Sample(matrix_.resample(mapping), transitions_.resample(selected, mapping));
    }

    friend std::ostream& operator<<(std::ostream &os, const Sample& o) { return o.print(os); }
    std::ostream& print(std::ostream &os) const {

        auto est_size = matrix_.num_features() * matrix_.num_states() * sizeof(FeatureMatrix::feature_value_t) /
                        (8 * 1024 * 1024);
        os << "Training sample [states: " << transitions_.num_states()
           << ", transitions: " << transitions_.num_transitions()
           << " (" << transitions_.num_marked_transitions() << " marked)"
           << ", goals: " << matrix_.num_goals()
           << ", features: " << matrix_.num_features()
           << ", est. size: " << est_size << " MB.]";
        return os;
    }
};



} // namespaces
