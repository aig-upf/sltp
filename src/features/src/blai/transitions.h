
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

class TransitionSample {
public:
    //! A transition is a pair of state IDs
    using transition_t = std::pair<unsigned, unsigned>;

    using transition_list_t = std::vector<transition_t>;
    using transition_set_t = std::unordered_set<transition_t, boost::hash<transition_t>>;

protected:
    const std::size_t num_states_;
    const std::size_t num_transitions_;
    const std::size_t num_marked_transitions_;
    // trdata_[s] contains the IDs of all neighbors of s in the state space
    std::vector<std::vector<unsigned>> trdata_;
    std::vector<bool> alive_states_;

    transition_set_t marked_transitions_;

public:
    TransitionSample(std::size_t num_states, std::size_t num_transitions, std::size_t num_marked_transitions)
            : num_states_(num_states),
              num_transitions_(num_transitions),
              num_marked_transitions_(num_marked_transitions),
              trdata_(num_states),
              alive_states_(num_states, false),
              marked_transitions_()
          {}

    ~TransitionSample() = default;
    TransitionSample(const TransitionSample&) = default;
    TransitionSample(TransitionSample&&) = default;

    std::size_t num_states() const { return num_states_; }
    std::size_t num_transitions() const { return num_transitions_; }
    std::size_t num_marked_transitions() const { return num_marked_transitions_; }

    const std::vector<unsigned>& successors(unsigned s) const {
        return trdata_.at(s);
    }

    const transition_set_t& marked_transitions() const {
        return marked_transitions_;
    }

    bool marked(const transition_t& p) const {
        return marked_transitions_.find(p) != marked_transitions_.end();
    }
    bool marked(unsigned src, unsigned dst) const {
        return marked(std::make_pair(src, dst));
    }

    bool is_alive(unsigned state) const {
        return alive_states_.at(state);
    }

    //! Print a representation of the object to the given stream.
    friend std::ostream& operator<<(std::ostream &os, const TransitionSample& o) { return o.print(os); }
    std::ostream& print(std::ostream &os) const {
        os << "Transition sample [states: " << num_states_ << ", transitions: " << num_transitions_;
        os << " (" << num_marked_transitions_ << " marked)]" << std::endl;
//        for (unsigned s = 0; s < num_states_; ++s) {
//            const auto& dsts = trdata_[s];
//            if (!dsts.empty()) os << "state " << s << ":";
//            for (auto dst:dsts) os << " " << dst;
//            os << std::endl;
//        }
        return os;
    }

    // readers
    void read(std::istream &is) {
        std::vector<transition_t> marked_transitions;
        marked_transitions.reserve(num_marked_transitions_);
        for( unsigned i = 0; i < num_marked_transitions_; ++i ) {
            unsigned src, dst;
            is >> src >> dst;
            assert(src < num_states_ && dst < num_states_);
            marked_transitions.emplace_back(src, dst);
        }

        // read number of states that have been expanded, for thich we'll have one state per line next
        unsigned num_records;
        is >> num_records;

        // read transitions, in format: source_id, num_successors, succ_1, succ_2, ...
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
                }
            }
        }

        // Validate that marked transitions are indeed transitions and store them
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

        // Store which states are alive (i.e. solvable, reachable, and not a goal)
        unsigned count, s;
        is >> count;
        assert(0 <= count && count <= num_states_);
        if(count > 0) {
            for(unsigned j = 0; j < count; ++j) {
                is >> s;
                assert(s < num_states_);
                alive_states_[s] = true;
            }
        }
    }

    static TransitionSample read_dump(std::istream &is, bool verbose) {
        unsigned num_states, num_transitions, num_marked_transitions;
        is >> num_states >> num_transitions >> num_marked_transitions;
        TransitionSample transitions(num_states, num_transitions, num_marked_transitions);
        transitions.read(is);
        if( verbose ) {
            std::cout << "TransitionSample::read_dump: #states=" << transitions.num_states()
                      << ", #transitions=" << transitions.num_transitions()
                      << ", #marked-transitions=" << transitions.marked_transitions_.size()
                      << std::endl;
        }
        return transitions;
    }

    //! Project the m states in selected, assumed to be a subset of [0, n], to the range
    //! [0, m], applying the given mapping
    TransitionSample resample(const std::unordered_set<unsigned>& selected,
            const std::unordered_map<unsigned, unsigned>& mapping) const {

        auto nstates = mapping.size();
        unsigned ntransitions = 0;

        std::vector<std::vector<unsigned>> trdata(nstates);
        transition_set_t marked_transitions;
        std::vector<bool> alive_states(nstates, false);

        for (unsigned s:selected) {
            unsigned mapped_s = mapping.at(s);
            assert(mapped_s < nstates);

            if (alive_states_.at(s)) alive_states.at(mapped_s) = true;

            for (unsigned sprime:successors(s)) {
                // mapping must contain all successors of the states in selected
                assert(mapping.find(sprime) != mapping.end());
                auto mapped_sprime = mapping.at(sprime);

                trdata.at(mapped_s).push_back(mapped_sprime);
                if (marked(s, sprime)) {
                    marked_transitions.emplace(mapped_s, mapped_sprime);
                }

                ++ntransitions;
            }
        }

        TransitionSample transitions(nstates, ntransitions, marked_transitions.size());
        transitions.trdata_ = std::move(trdata);
        transitions.marked_transitions_ = std::move(marked_transitions);
        return transitions;
    }
};

} // namespaces
