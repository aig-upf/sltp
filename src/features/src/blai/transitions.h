
#pragma once

#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>


namespace Sample {

class Transitions {
public:
    //! A transition is a pair of state IDs
    using transition_t = std::pair<unsigned, unsigned>;

    using transition_list_t = std::vector<transition_t>;
    using transition_set_t = std::unordered_set<transition_t, boost::hash<transition_t>>;

protected:
    const unsigned num_states_;
    const unsigned num_transitions_;
    unsigned num_marked_transitions_;
    // trdata_[s] contains the IDs of all neighbors of s in the state space
    std::vector<std::vector<unsigned>> trdata_;
    transition_list_t all_transitions_;

    transition_set_t marked_transitions_;

public:
    Transitions(unsigned num_states, unsigned num_transitions, unsigned num_marked_transitions)
            : num_states_(num_states),
              num_transitions_(num_transitions),
              num_marked_transitions_(num_marked_transitions),
              trdata_(num_states)
          {}
    ~Transitions() = default;

    unsigned num_states() const {
        return num_states_;
    }
    unsigned num_transitions() const {
        return num_transitions_;
    }
    std::size_t num_transitions(unsigned s) const {
        return trdata_.at(s).size();
    }

    const std::vector<unsigned>& successors(unsigned s) const {
        return trdata_.at(s);
    }

    const transition_set_t& marked_transitions() const {
        return marked_transitions_;
    }

    const transition_list_t& all_transitions() const {
        return all_transitions_;
    }

    bool marked(const transition_t& p) const {
        return marked_transitions_.find(p) != marked_transitions_.end();
    }
    bool marked(unsigned src, unsigned dst) const {
        return marked(std::make_pair(src, dst));
    }

    void print(std::ostream &os) const {
        os << "Transition set: #states=" << num_states_ << ", #transitions=" << num_transitions_ << std::endl;
        for (unsigned s = 0; s < num_states_; ++s) {
            const auto& dsts = trdata_[s];
            if (!dsts.empty()) os << "state " << s << ":";
            for (auto dst:dsts) os << " " << dst;
            os << std::endl;
        }
    }

    // readers
    void read(std::istream &is) {
        // read marked transitions (this is somewhat redundant, but it was added afterwards)
        unsigned num_marked_transitions;
        is >> num_marked_transitions;
        assert(num_marked_transitions == num_marked_transitions_);

        std::vector<transition_t> marked_transitions;
        marked_transitions.reserve(num_marked_transitions);
        for( unsigned i = 0; i < num_marked_transitions_; ++i ) {
            unsigned src, dst;
            is >> src >> dst;
            assert(src < num_states_ && dst < num_states_);
            marked_transitions.emplace_back(src, dst);
        }

        // read number of records in the rest of the file
        unsigned num_records;
        is >> num_records;

        // read transitions
        for( unsigned i = 0; i < num_records; ++i ) {
            unsigned src, count, dst;
            is >> src >> count;
            assert(src < num_states_ && 0 <= count);
            if( count > 0 ) {
                std::vector<bool> seen(num_states_, false);
                trdata_[src].reserve(count);
                for( unsigned j = 0; j < count; ++j ) {
                    is >> dst;
                    assert(dst < num_states_);
                    if (seen.at(dst)) throw std::runtime_error("Duplicate transition");
                    trdata_[src].push_back(dst);
                    seen[dst] = true;
                    all_transitions_.emplace_back(src, dst);
                }
            }
        }

        // store valid marked transitions
        for (const auto &marked : marked_transitions) {
            unsigned src = marked.first;
            unsigned dst = marked.second;

            // Check that the marked transition is indeed a transition
            bool valid = false;
            for (unsigned t:trdata_[src]) {
                if (dst == t) {
                    valid = true;
                    break;
                }
            }

            if (!valid) {
                throw std::runtime_error("Invalid marked transition");
            }
            marked_transitions_.emplace(src, dst);
        }
    }
    static Transitions read_dump(std::istream &is, bool verbose) {
        unsigned num_states, num_transitions, num_marked_transitions;
        is >> num_states >> num_transitions >> num_marked_transitions;
        Transitions transitions(num_states, num_transitions, num_marked_transitions);
        transitions.read(is);
        if( verbose ) {
            std::cout << "Transitions::read_dump: #states=" << transitions.num_states()
                      << ", #transitions=" << transitions.num_transitions()
                      << ", #marked-transitions=" << transitions.marked_transitions_.size()
                      << std::endl;
        }
        return transitions;
    }
};

} // namespaces
