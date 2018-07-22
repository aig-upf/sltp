/*
This program computes the distance to goal (i.e. the cost of the least-cost path to the goal)
of every state from which the goal can be reached.
It does this by executing Dijkstra's algorithm backwards from the goal.
It prints on stdout each state and its distance (distance first, then the state) and, if a filename is
provided as a command line argument, it prints the state_map it builds to that file.

Copyright (C) 2013 by the PSVN Research Group, University of Alberta
*/

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "priority_queue.hpp"

using namespace std;

class Feature {
  protected:
    mutable int value_;
  public:
    Feature() : value_(0) { }
    virtual ~Feature() { }
    virtual string name() const = 0;
    virtual void calculate(int num_blocks, const state_t &state) const = 0;
    int value() const { return value_; }
};

class BooleanFeature : public Feature {
  public:
    BooleanFeature() { }
    virtual ~BooleanFeature() { }
};

class NumericalFeature : public Feature {
  public:
    NumericalFeature() { }
    virtual ~NumericalFeature() { }
    int num_above(int block, int num_blocks, const state_t &state) const {
        int count = 0;
        if( state.vars[block] == 0 ) { // not clear(block)
            for( int b = 1; b <= num_blocks; ++b ) {
                if( state.vars[b] == 1 ) {
                    int t = b;
                    while( (t != 0) && (t != block) ) {
                        ++count;
                        t = state.vars[num_blocks + t];
                    }
                    if( t == block ) break;
                }
                count = 0;
            }
            assert(count > 0);
        }
        return count;
    }
    int num_below(int block, int num_blocks, const state_t &state) const {
        int count = 0;
        int b = state.vars[num_blocks + block];
        while( b != 0 ) {
            ++count;
            b = state.vars[num_blocks + b];
        }
        return count;
    }
    int num_other(int block, int num_blocks, const state_t &state) const {
        int nabove = num_above(block, num_blocks, state);
        int nbelow = num_below(block, num_blocks, state);
        assert(1 + nbelow + nabove <= num_blocks);
        return num_blocks - nbelow - nabove - 1;
    }
};


class DummyFeature : public BooleanFeature {
  public:
    DummyFeature() { }
    virtual ~DummyFeature() { }
    virtual string name() const {
        return string("dummy()");
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = 0;
    }
};

class HoldingBlock : public BooleanFeature {
  protected:
    int block_;
  public:
    HoldingBlock(int block) : block_(block) { }
    virtual ~HoldingBlock() { }
    virtual string name() const {
        return block_ == 0 ? string("empty()") : string("hold(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = state.vars[0] == block_;
    }
};

class HoldingOther : public BooleanFeature {
  protected:
    int block_;
  public:
    HoldingOther(int block) : block_(block) { }
    virtual ~HoldingOther() { }
    virtual string name() const {
        return string("hold-other(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = (state.vars[0] != 0) && (state.vars[0] != block_);
    }
};

class ClearBlock : public BooleanFeature {
  protected:
    int block_;
  public:
    ClearBlock(int block) : block_(block) { }
    virtual ~ClearBlock() { }
    virtual string name() const {
        return string("clear(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = state.vars[block_];
    }
};

class OnTable : public BooleanFeature {
  protected:
    int block_;
  public:
    OnTable(int block) : block_(block) { }
    virtual ~OnTable() { }
    virtual string name() const {
        return string("ontable(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = state.vars[num_blocks + block_] == 0;
    }
};

class SomeBelow : public BooleanFeature { // this is opposite of ontable(x)
  protected:
    int block_;
  public:
    SomeBelow(int block) : block_(block) { }
    virtual ~SomeBelow() { }
    virtual string name() const {
        return string("some-below(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = state.vars[num_blocks + block_] != 0;
    }
};

class Above : public BooleanFeature {
  protected:
    int block1_;
    int block2_;
  public:
    Above(int block1, int block2) : block1_(block1), block2_(block2) { }
    virtual ~Above() { }
    virtual string name() const {
        return string("above(") + domain_0[block1_] + "," + domain_0[block2_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        int block = block1_;
        while( (block != 0) && (block != block2_) )
            block = state.vars[num_blocks + block];
        value_ = block == block2_;
    }
};

class Below : public BooleanFeature {
  protected:
    int block1_;
    int block2_;
  public:
    Below(int block1, int block2) : block1_(block1), block2_(block2) { }
    virtual ~Below() { }
    virtual string name() const {
        return string("below(") + domain_0[block1_] + "," + domain_0[block2_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        int block = block2_;
        while( (block != 0) && (block != block1_) )
            block = state.vars[num_blocks + block];
        value_ = block == block1_;
    }
};

class On : public BooleanFeature {
  protected:
    int block1_;
    int block2_;
  public:
    On(int block1, int block2) : block1_(block1), block2_(block2) { }
    virtual ~On() { }
    virtual string name() const {
        return string("on(") + domain_0[block1_] + "," + domain_0[block2_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = state.vars[num_blocks + block1_] == block2_;
    }
};

class NumBlocksAbove : public NumericalFeature {
  protected:
    int block_;
  public:
    NumBlocksAbove(int block) : block_(block) { }
    virtual ~NumBlocksAbove() { }
    virtual string name() const {
        return string("nabove(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = num_above(block_, num_blocks, state);
    }
};

class NumBlocksBelow : public NumericalFeature {
  protected:
    int block_;
  public:
    NumBlocksBelow(int block) : block_(block) { }
    virtual ~NumBlocksBelow() { }
    virtual string name() const {
        return string("nbelow(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        value_ = num_below(block_, num_blocks, state);
    }
};

class NumBlocksOther : public NumericalFeature {
  protected:
    int block_;
  public:
    NumBlocksOther(int block) : block_(block) { }
    virtual ~NumBlocksOther() { }
    virtual string name() const {
        return string("nother(") + domain_0[block_] + ")";
    }
    virtual void calculate(int num_blocks, const state_t &state) const {
        if( state.vars[0] != block_ ) {
            value_ = num_other(block_, num_blocks, state);
            value_ -= state.vars[0] != 0 ? 1 : 0; // subtract 1 if holding some block
        } else {
            value_ = num_blocks - 1;
        }
    }
};

class FeaturesAndTransitions {
  protected:
    struct StateComparator {
      int operator()(const state_t &s1, const state_t &s2) const {
          return compare_states(&s1, &s2) < 0;
      }
    };

  protected:
    int num_blocks_;
    int num_transitions_;
    int last_numerical_feature_;
    int first_boolean_feature_;
    vector<const Feature*> features_;
    map<state_t, int, StateComparator> states_;
    map<int, set<pair<int, int> > > transitions_;

  public:
    FeaturesAndTransitions(int num_blocks)
      : num_blocks_(num_blocks), num_transitions_(0), last_numerical_feature_(0) { }
    ~FeaturesAndTransitions() {
        for( size_t i = 0; i < features_.size(); ++i )
            delete features_[i];
    }

    int num_states() const {
        return states_.size();
    }
    int num_transitions() const {
        return num_transitions_;
    }
    int num_features() const {
        return features_.size();
    }
    int last_numerical_feature() const {
        return last_numerical_feature_;
    }
    int first_boolean_feature() const {
        return first_boolean_feature_;
    }

    void create_features() {
        // numerical features
        for( int block = 1; block <= num_blocks_; ++block ) {
            features_.push_back(new NumBlocksAbove(block));
            features_.push_back(new NumBlocksBelow(block));
            features_.push_back(new NumBlocksOther(block));
        }
        last_numerical_feature_ = features_.size();

        // calculate index first boolean feature
        if( last_numerical_feature_ == 0 ) {
            first_boolean_feature_ = 0;
        } else {
            first_boolean_feature_ = 1;
            while( first_boolean_feature_ < last_numerical_feature_ )
                first_boolean_feature_ = first_boolean_feature_ << 1;
        }

        // insert dummy features until first boolean feature
        while( features_.size() < first_boolean_feature_ )
            features_.push_back(new DummyFeature());

        // boolean features
        features_.push_back(new HoldingBlock(0)); // holding nothing
        for( int block = 1; block <= num_blocks_; ++block ) {
            features_.push_back(new HoldingBlock(block));
            features_.push_back(new HoldingOther(block));
            features_.push_back(new ClearBlock(block));
            features_.push_back(new OnTable(block));
            features_.push_back(new SomeBelow(block));
        }
        for( int block1 = 1; block1 <= num_blocks_; ++block1 ) {
            for( int block2 = 1; block2 <= num_blocks_; ++block2 ) {
                if( block1 != block2 ) {
                    features_.push_back(new Above(block1, block2));
                    features_.push_back(new Below(block1, block2));
                    features_.push_back(new On(block1, block2));
                }
            }
        }
    }

    void print_features(ostream &os) const {
        os << features_.size();
        for( size_t i = 0; i < features_.size(); ++i )
            os << " " << features_[i]->name();
        os << endl;
    }

    void calculate_features(const state_t &state) const {
        for( size_t i = 0; i < features_.size(); ++i )
            features_[i]->calculate(num_blocks_, state);
    }

    void dump_feature_valuation(ostream &os, const state_t &state, int state_index, bool boolean, bool names) const {
        //print_state(stdout, &state); printf(" ");
        int count = 0;
        for( size_t i = 0; i < features_.size(); ++i )
            count += features_[i]->value() > 0 ? 1 : 0;
        os << state_index << " " << count;

        for( size_t i = 0; i < features_.size(); ++i ) {
            if( features_[i]->value() > 0 ) {
                if( names )
                    os << " " << features_[i]->name();
                else
                    os << " " << i;
                if( !boolean )
                    os << ":" << features_[i]->value();
            }
        }
        os << endl;
    }
    void dump_features_for_states(ostream &os, bool boolean, bool names) const {
        for( map<state_t, int>::const_iterator it = states_.begin(); it != states_.end(); ++it ) {
            const state_t &state = it->first;
            int state_index = it->second;
            calculate_features(state);
            dump_feature_valuation(os, state, state_index, boolean, names);
        }
    }
    void dump_features(ostream &os, bool boolean = true, bool names = false) const {
        os << num_states() << " "
           << num_features() << " "
           << last_numerical_feature() << " "
           << first_boolean_feature()
           << endl;
        print_features(os);
        dump_features_for_states(os, boolean, names);
    }

    void add_transition(int src, int label, int dst) {
        pair<int, int> p(dst, label);
        map<int, set<pair<int, int> > >::const_iterator it = transitions_.find(src);
        if( (it == transitions_.end()) || (it->second.find(p) == it->second.end()) ) {
            ++num_transitions_;
            //cout << "new transition: " << src << " " << dst << endl;
        }
        transitions_[src].insert(p);
    }
    void add_transition(const state_t &src, int label, const state_t &dst) {
        pair<map<state_t, int>::iterator, bool> p_dst = states_.insert(make_pair(dst, states_.size()));
        pair<map<state_t, int>::iterator, bool> p_src = states_.insert(make_pair(src, states_.size()));
        add_transition(p_src.first->second, label, p_dst.first->second);
    }

    void dump_transitions(ostream &os, bool labels = false) const {
        os << num_states() << " " << num_transitions() << endl;
        for( int src = 0; src < num_states(); ++src ) {
            map<int, set<pair<int, int> > >::const_iterator it = transitions_.find(src);
            if( it == transitions_.end() ) {
                os << src << " 0";
            } else {
                assert(!it->second.empty());
                os << src << " " << it->second.size();
                for( set<pair<int, int> >::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt ) {
                    int dst = jt->first;
                    int label = jt->second;
                    if( labels ) os << " " << get_bwd_rule_label(label);
                    os << " " << dst;
                }
            }
            os << endl;
        }
    }

    void dump_goal_features(ostream &os) const {
        vector<int> goal_features;
        for( size_t i = 0; i < features_.size(); ++i ) {
            if( features_[i]->name() == "nabove(B1)" )
                goal_features.push_back(i);
            else if( features_[i]->name() == "hold(B1)" )
                goal_features.push_back(i);
        }
        os << goal_features.size();
        for( size_t i = 0; i < goal_features.size(); ++i )
            os << " " << goal_features[i];
        os << endl;
    }
};

int main(int argc, char **argv) {
    // usage
    if( argc < 3 ) {
        cout << "Usage: " << argv[0] << " <num-blocks> <boolean>" << endl;
        exit(0);
    }

    // features
    int num_blocks = atoi(argv[1]);
    bool boolean = atoi(argv[2]) == 1;
    cout << "parameters: #blocks=" << num_blocks << ", boolean=" << boolean << endl;

    FeaturesAndTransitions FT(num_blocks);
    FT.create_features();

    PriorityQueue<state_t> open; // used for the states we have generated but not yet expanded (the OPEN list)
    state_map_t *state_map = new_state_map(); // contains the cost-to-goal for all states that have been generated

    // add goal states
    int d = 0;
    state_t state;
    first_goal_state(&state, &d);
    do {
        state_map_add(state_map, &state, 0);
        open.Add(0, 0, state);
    } while( next_goal_state(&state, &d) );

    while( !open.Empty() ) {
        // get current distance from goal
        d = open.CurrentPriority();

        // get state
        state_t state = open.Top();
        open.Pop();
        
        // check if we already expanded this state.
	// (entries on the open list are not deleted if a cheaper path to a state is found)
        const int *best_dist = state_map_get(state_map, &state);
        assert(best_dist != nullptr);
        if( *best_dist < d ) continue;
        
        // print the distance then the state
        //printf("%d  ",d);
        //print_state(stdout, &state);
        //printf(" \n");

        // look at all predecessors of the state
        ruleid_iterator_t iter;
        init_bwd_iter(&iter, &state);
        for( int ruleid = next_ruleid(&iter); ruleid >= 0; ruleid = next_ruleid(&iter) ) {
            state_t child;   // NOTE: "child" will be a predecessor of state, not a successor
            apply_bwd_rule(ruleid, &state, &child);
            const int child_d = d + get_bwd_rule_cost(ruleid);

            // check if either this child has not been seen yet or if
            // there is a new cheaper way to get to this child.
            const int *old_child_d = state_map_get(state_map, &child);
            if( (old_child_d == nullptr) || (*old_child_d > child_d) ) {
                // add to open with the new distance
                state_map_add(state_map, &child, child_d);
                open.Add(child_d, child_d, child);
            }

            // add transition
            FT.add_transition(child, ruleid, state);
        }
    }

    string prefix = string("blocks") + (num_blocks < 10 ? "0" : "") + to_string(num_blocks);
    string suffix = string(".dat") + (boolean ? "" : "2");
    ofstream osm((prefix + "_matrix" + suffix).c_str());
    FT.dump_features(osm, boolean);
    osm.close();
    ofstream ost((prefix + "_transitions" + suffix).c_str());
    FT.dump_transitions(ost);
    ost.close();
    ofstream osg((prefix + "_goal_features" + suffix).c_str());
    FT.dump_goal_features(osg);
    osg.close();
   
    /* 
    // write the state map to a file
    if( argc >= 2 ) {
        FILE *file = fopen(argv[1], "w");
        if( file == nullptr ) {
            fprintf(stderr, "could not open %s for writing\n", argv[1]);
            exit(-1);
        }
        write_state_map(file, state_map);
        fclose(file);
    }
    */
    
    return 0;
}
