#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <theory.h>
#include "sample.h"

// TODO: * action must be applicable at some state
//       * further ordering of actions to reduce more symmetries
//       * separate goal from non-goal states (IMPORTANT)

using namespace std;

string filename(const string &prefix, int K, int N, const string &suffix, bool opt = true) {
    string fn(prefix);
    if( opt ) {
        fn += string("_") + to_string(K);
        fn += string("_") + to_string(N);
    }
    fn += suffix;
    return fn;
}

template<typename T>
class Items : public set<T> {
    int num_;
  public:
    Items(int num = 0) : num_(num) { }
    ~Items() { }
    int num() const {
        return num_;
    }
    void dump(ostream &os) const {
        os << set<T>::size();
        for( typename set<T>::const_iterator it = set<T>::begin(); it != set<T>::end(); ++it )
            os << " " << *it;
        os << endl;
    }
    void read(istream &is) {
        for( size_t i = 0; i < num_; ++i ) {
            T item;
            is >> item;
            set<T>::insert(item);
        }
    }
    static const Items<T>* read_dump(istream &is) {
        int num;
        is >> num;
        Items<T> *items = new Items<T>(num);
        items->read(is);
        cout << "Items::read_dump: #items=" << items->num() << endl;
        return items;
    }
};

class Theory : public SAT::Theory {
  public:
    struct AbstractAction {
        string name_;
        int selected_precondition_;
        int precondition_;
        int selected_effect_;
        int effect_;

        static map<string, int> *feature_map_;

        AbstractAction() : selected_precondition_(0), precondition_(0), selected_effect_(0), effect_(0) { }
        AbstractAction(const string &name) : name_(name), selected_precondition_(0), precondition_(0), selected_effect_(0), effect_(0) { }

        void add_precondition(const string &feature, int prec) {
            if( (feature_map_ == nullptr) || (feature_map_->find(feature) == feature_map_->end()) ) {
                cout << "error: add_precondition: feature '" << feature << "' not found in feature map; did you forget --add-features to add features in forced actions?" << endl;
                exit(0);
            }
            assert((prec == 0) || (prec == 1));
            int k = feature_map_->find(feature)->second;
            selected_precondition_ = selected_precondition_ | (1 << k);
            precondition_ = precondition_ & ~(1 << k);
            precondition_ += prec << k;
        }
        void add_effect(const string &feature, int effect) {
            if( (feature_map_ == nullptr) || (feature_map_->find(feature) == feature_map_->end()) ) {
                cout << "error: add_precondition: feature '" << feature << "' not found in feature map; did you forget --add-features to add features in forced actions?" << endl;
                exit(0);
            }
            assert((effect == 0) || (effect == 1));
            int k = feature_map_->find(feature)->second;
            selected_effect_ = selected_effect_ | (1 << k);
            effect_ = effect_ & ~(1 << k);
            effect_ += effect << k;
        }
        bool operator<(const AbstractAction &a) const {
            return
              (selected_precondition_ < a.selected_precondition_) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ < a.precondition_)) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ == a.precondition_) && (selected_effect_ < a.selected_effect_)) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ == a.precondition_) && (selected_effect_ == a.selected_effect_) && (effect_ < a.effect_));
        }
        void print(ostream &os, const Theory &T, const vector<bool> &model, bool qnp_format) const {
            if( qnp_format ) {
                os << name_ << endl;

                int num_preconditions = 0;
                for( int k = 0; k < T.K_; ++k )
                    num_preconditions += selected_precondition_ & (1 << k) ? 1 : 0;
                os << num_preconditions;
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_precondition_ & (1 << k) ) {
                        os << " " << T.condition_on_selected_feature(model, k, (precondition_ & (1 << k)) != 0, qnp_format);
                    }
                }
                os << endl;

                int num_effects = 0;
                for( int k = 0; k < T.K_; ++k )
                    num_effects += selected_effect_ & (1 << k) ? 1 : 0;
                os << num_effects;
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_effect_ & (1 << k) ) {
                        os << " " << T.effect_on_selected_feature(model, k, (effect_ & (1 << k)) != 0, qnp_format);
                    }
                }
                os << endl;
            } else {
                os << "Prec={";
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_precondition_ & (1 << k) )
                        os << " " << T.condition_on_selected_feature(model, k, (precondition_ & (1 << k)) != 0, qnp_format) << ",";
                }
                os << " }; Effect={";
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_effect_ & (1 << k) )
                        os << " " << T.effect_on_selected_feature(model, k, (effect_ & (1 << k)) != 0) << ",";
                }
                os << " }" << endl;
            }
        }
        void dump(ostream &os) const {
            os << SAT::bit_string(selected_precondition_)
               << " " << SAT::bit_string(precondition_)
               << " " << SAT::bit_string(selected_effect_)
               << " " << SAT::bit_string(effect_)
               << endl;
        }
    };

    struct Options {
        string prefix_;
        bool ordering_;
        bool decode_;
        bool explicit_goals_;
        bool add_actions_;
        bool add_features_;
        bool only_soundness_;
        bool log_feature_variables_;
        Options()
          : ordering_(false),
            decode_(false),
            explicit_goals_(false),
            add_actions_(false),
            add_features_(false),
            only_soundness_(false),
            log_feature_variables_(true) {
        }
    };

  protected:
    const Sample::Sample &S_;
    const int K_;           // number of features in abstraction
    const int N_;           // number of actions in abstraction
    const Options options_;

    // forced elements (added by "hand")
    Items<AbstractAction> forced_actions_;
    Items<string> forced_features_;
    Items<int> forced_goals_;

    const int num_states_;
    const int num_features_;
    int log_num_features_;
    int log_last_numerical_feature_;

  public:
    Theory(const Sample::Sample &S, int K, int N, const Options &options)
      : SAT::Theory(options.decode_),
        S_(S), K_(K), N_(N),
        options_(options),
        num_states_(S_.M().num_states()),
        num_features_(S_.M().num_features()),
        log_num_features_(0),
        log_last_numerical_feature_(0) {

        while( (1 << log_num_features_) < num_features_ )
            ++log_num_features_;
        while( (1 << log_last_numerical_feature_) < S_.M().last_numerical_feature() )
            ++log_last_numerical_feature_;

        cout << "Theory: parameters:"
             << " K=" << K_
             << ", N=" << N_
             << ", log2(#features)=" << log_num_features_
             << ", log2(last_numerical_feature)=" << log_last_numerical_feature_
             << endl;

        if( options_.add_features_ )
            read_features(options_.prefix_);
        if( options_.explicit_goals_ )
            read_goals(options_.prefix_);
        if( options_.add_actions_ )
            read_actions(options_.prefix_);

        build_theory();
    }
    virtual ~Theory() { }

    virtual void build_rest() {
        add_implications_for_features();
        add_implications_for_goals();
        add_implications_for_actions();

        cout << "------ soundness + completeness ------" << endl;
        enforce_soundness();
        if( !options_.only_soundness_ ) enforce_completeness();
        cout << "--------------------------------------" << endl;
    }

    int index_for_selected_feature(const vector<bool> &model, int k) const {
        int index = 0;
        if( options_.log_feature_variables_ ) {
            for( int t = 0; t < log_num_features_; ++t ) {
                assert(F(k, t) < model.size());
                index += model[F(k, t)] ? (1 << t) : 0;
            }
        } else {
            index = -1;
            for( int j = 0; j < num_features_; ++j ) {
                if( model[Fn(k, j)] ) {
                    index = j;
                    break;
                }
            }
            assert(index != -1);
        }
        return index;
    }
    string condition_on_selected_feature(const vector<bool> &model, int k, bool prec, bool qnp_format = false) const {
        // get feature index from model
        int f = index_for_selected_feature(model, k);

        // generate string
        string feature;
        if( qnp_format ) {
            feature += S_.M().feature(f);
            feature += prec ? " 1" : " 0";
        } else if( S_.M().numerical(f) ) {
            feature += to_string(f) + "." + S_.M().feature(f);
            feature += prec ? ">0" : "=0";
        } else {
            feature += prec ? "" : "-";
            feature += to_string(f) + "." + S_.M().feature(f);
        }
        return feature;
    }
    string effect_on_selected_feature(const vector<bool> &model, int k, bool eff, bool qnp_format = false) const {
        // get feature index from model
        int f = index_for_selected_feature(model, k);

        // generate string
        string feature;
        if( qnp_format ) {
            feature += S_.M().feature(f);
            feature += eff ? " 1" : " 0";
        } else if( S_.M().numerical(f) ) {
            feature += eff ? "INC(" : "DEC(";
            feature += to_string(f) + "." + S_.M().feature(f) + ")";
        } else {
            feature += eff ? "" : "-";
            feature += to_string(f) + "." + S_.M().feature(f);
        }
        return feature;
    }

  protected:
    bool do_ordering() const {
        bool disable = (forced_features_.size() == K_) || ((forced_features_.size() < K_) && !forced_actions_.empty());
        return options_.ordering_ && !disable;
    }

    virtual void build_variables() {
        // A variables: 4 x K x N
        var_offsets_.push_back(make_pair(0, "A"));
        for( int t = 1; t < 5; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    int index = variables_.size();
                    string name = string("A(") + to_string(t) + "," + to_string(j) + "," + to_string(k) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == A(t, j, k));
                }
            }
        }
        cout << "A: #variables=" << variables_.size() - var_offsets_.back().first << " (4KN), offset=" << var_offsets_.back().first << endl;

        if( options_.log_feature_variables_ ) {
            // F variables: K x log2(F) (logarithmic encoding)
            var_offsets_.push_back(make_pair(variables_.size(), "F"));
            for( int k = 0; k < K_; ++k ) {
                for( int t = 0; t < log_num_features_; ++t ) {
                    int index = variables_.size();
                        string name = string("F(") + to_string(k) + "," + to_string(t) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == F(k, t));
                }
            }
            cout << "F: #variables=" << variables_.size() - var_offsets_.back().first << " (KL) where L=log2(F), offset=" << var_offsets_.back().first << endl;
        } else {
            // Fn variables: K x F (naive encoding)
            var_offsets_.push_back(make_pair(variables_.size(), "Fn"));
            for( int k = 0; k < K_; ++k ) {
                for( int j = 0; j < num_features_; ++j ) {
                    int index = variables_.size();
                    string name = string("Fn(") + to_string(k) + "," + to_string(j) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == Fn(k, j));
                }
            }
            cout << "Fn: #variables=" << variables_.size() - var_offsets_.back().first << " (KF), offset=" << var_offsets_.back().first << endl;
        }

        // n variables: K
        var_offsets_.push_back(make_pair(variables_.size(), "n"));
        for( int k = 0; k < K_; ++k ) {
            int index = variables_.size();
            string name = string("n(") + to_string(k) + ")";
            variables_.push_back(new SAT::Var(index, name));
            assert(index == n(k));
        }
        cout << "n: #variables=" << variables_.size() - var_offsets_.back().first << " (K), offset=" << var_offsets_.back().first << endl;

        // phi and phi* variables: 2 x S x K
        var_offsets_.push_back(make_pair(variables_.size(), "phi"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                int index = variables_.size();
                string name = string("phi(") + to_string(k) + "," + to_string(i) + ")";
                variables_.push_back(new SAT::Var(index, name));
                assert(index == phi(k, i, false));
            }
        }
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                int index = variables_.size();
                string name = string("phi*(") + to_string(k) + "," + to_string(i) + ")";
                variables_.push_back(new SAT::Var(index, name));
                assert(index == phi(k, i, true));
            }
        }
        cout << "phi: #variables=" << variables_.size() - var_offsets_.back().first << " (2SK), offset=" << var_offsets_.back().first << endl;

        // B variables: 2 x S x K x N
        var_offsets_.push_back(make_pair(variables_.size(), "B"));
        for( int t = 1; t < 3; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int i = 0; i < num_states_; ++i ) {
                    for( int k = 0; k < K_; ++k ) {
                        int index = variables_.size();
                        string name = string("B(") + to_string(t) + "," + to_string(j) + "," + to_string(i) + "," + to_string(k) + ")";
                        variables_.push_back(new SAT::Var(index, name));
                        assert(index == B(t, j, i, k));
                    }
                }
            }
        }
        cout << "B: #variables=" << variables_.size() - var_offsets_.back().first << " (2SKN), offset=" << var_offsets_.back().first << endl;

        // RES variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "RES"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("RES(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == RES(i, l, j));
                }
            }
        }
        cout << "RES: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        // INC variables: T x K
        var_offsets_.push_back(make_pair(variables_.size(), "INC"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("INC(") + to_string(k) + "," + to_string(i) + "," + to_string(l) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == INC(k, i, l));
                }
            }
        }
        cout << "INC: #variables=" << variables_.size() - var_offsets_.back().first << " (TK), offset=" << var_offsets_.back().first << endl;

        // DEC variables: T x K
        var_offsets_.push_back(make_pair(variables_.size(), "DEC"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("DEC(") + to_string(k) + "," + to_string(i) + "," + to_string(l) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == DEC(k, i, l));
                }
            }
        }
        cout << "DEC: #variables=" << variables_.size() - var_offsets_.back().first << " (TK), offset=" << var_offsets_.back().first << endl;

        // MATCH variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "MATCH"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("MATCH(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == MATCH(i, l, j));
                }
            }
        }
        cout << "MATCH: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        // APPMATCH variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "APPMATCH"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("APPMATCH(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == APPMATCH(i, l, j));
                }
            }
        }
        cout << "APPMATCH: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        if( do_ordering() ) {
            if( options_.log_feature_variables_ ) {
                // df variables:
                var_offsets_.push_back(make_pair(variables_.size(), "df"));
                for( int j = 0; j + 1 < K_; ++j ) {
                    for( int k = 0; k < log_num_features_; ++k ) {
                        int index = variables_.size();
                        string name = string("df(") + to_string(j) + "," + to_string(k) + ")";
                        variables_.push_back(new SAT::Var(index, name));
                        assert(index == df(j, k));
                    }
                }
                cout << "df: #variables=" << variables_.size() - var_offsets_.back().first << " ((K-1) x log2(F)), offset=" << var_offsets_.back().first << endl;
            } else {
                // there is no need for extra variables when log encoding is disabled
            }
        }

        add_variables_for_features();
        add_variables_for_goals();
    }

    // references to variables
    // A(1, j, k) iff j-th action *affects* k-th feature
    // A(2, j, k) iff type of effect of j-th action on k-th feature (0=false/decrement, 1=true/increment)
    // A(3, j, k) iff k-th feature is *precondition* of j-th action
    // A(4, j, k) iff type of precondition of k-th feature in j-th action (0=false/=0, 1=true/>0)
    int A(int t, int j, int k) const { // t in {1, ... ,4}, j in {0, ..., N-1}, k in {0, ... ,K-1}
        assert((0 < t) && (t < 5));
        assert((0 <= j) && (j < N_));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * (j + N_ * (t - 1));
        assert(var_offsets_[0].second == "A");
        return var_offsets_[0].first + index;
    }

    // F(k, t) iff k-th feature is feature whose index has t-th bit set to 1 (logarithmic encoding)
    int F(int k, int t) const { // k in {0, ..., K-1}, t in {0, ..., log2(F)-1}
        assert(options_.log_feature_variables_);
        assert((0 <= k) && (k < K_));
        assert((0 <= t) && (t < log_num_features_));
        int index = t + log_num_features_ * k;
        assert(var_offsets_[1].second == "F");
        return var_offsets_[1].first + index;
    }
    // Fn(k, j) iff k-th feature is feature whose index is j (naive encoding)
    int Fn(int k, int j) const { // k in {0, ..., K-1}, j in {0, ..., F-1}
        assert(!options_.log_feature_variables_);
        assert((0 <= k) && (k < K_));
        assert((0 <= j) && (j < num_features_));
        int index = j + num_features_ * k;
        assert(var_offsets_[1].second == "Fn");
        return var_offsets_[1].first + index;
    }

    int n(int k) const { // k in {0, ..., K-1}
        assert((0 <= k) && (k < K_));
        int index = k;
        assert(var_offsets_[2].second == "n");
        return var_offsets_[2].first + index;
    }

    int df(int k, int j) const { // k in {0, ..., K-2}, j in {0, ..., log2(F)-1}
        assert((0 <= k) && (k + 1 < K_));
        assert((0 <= j) && (j < log_num_features_));
        int index = j + log_num_features_ * k;
        assert(var_offsets_[10].second == "df");
        return var_offsets_[10].first + index;
    }
    int MappedBy(int j, int k) const { // j in {0, ..., #forced-features-1}, k in {0, ..., K-1}
        assert((0 <= j) && (j < int(forced_features_.size())));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * j;
        assert(var_offsets_[11].second == "MappedBy");
        return var_offsets_[11].first + index;
    }

    // phi(k, i, F) iff k-th feature *holds* in state s_i in sample
    // phi(k, i, T) iff k-th feature *does not hold* in state s_i in sample
    int phi(int k, int i, bool star) const { // k in {0, ..., K-1}, i in {0, ..., states-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        int index = i + num_states_ * (k + K_ * int(star));
        assert(var_offsets_[3].second == "phi");
        return var_offsets_[3].first + index;
    }

    // used for soundness
    int B(int t, int j, int i, int k) const { // t in {1,...,2}, j in {0, ..., N-1}, i in {0, ..., states-1}, k in {0, ..., K-1}
        assert((1 <= t) && (t < 3));
        assert((0 <= j) && (j < N_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * (i + num_states_ * (j + N_ * (t - 1)));
        assert(var_offsets_[4].second == "B");
        return var_offsets_[4].first + index;
    }

    // used for soundness/completeness
    int RES(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (S_.T().offset(i) + l) + S_.T().num_transitions() * j;
        assert(var_offsets_[5].second == "RES");
        return var_offsets_[5].first + index;
    }

    // INC(k, i, l) iff F^k increases in transition (s_i, s_l) in sample
    int INC(int k, int i, int l) const { // k in {0, ..., K}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        int index = (S_.T().offset(i) + l) + S_.T().num_transitions() * k;
        assert(var_offsets_[6].second == "INC");
        return var_offsets_[6].first + index;
    }

    // DEC(k, i, l) iff F^k decreases in transition (s_i, s_l) in sample
    int DEC(int k, int i, int l) const { // k in {0, ..., K}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        int index = (S_.T().offset(i) + l) + S_.T().num_transitions() * k;
        assert(var_offsets_[7].second == "DEC");
        return var_offsets_[7].first + index;
    }

    // used for soundness/completeness
    int MATCH(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (S_.T().offset(i) + l) + S_.T().num_transitions() * j;
        assert(var_offsets_[8].second == "MATCH");
        return var_offsets_[8].first + index;
    }

    // used only for completeness
    int APPMATCH(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (S_.T().offset(i) + l) + S_.T().num_transitions() * j;
        assert(var_offsets_[9].second == "APPMATCH");
        return var_offsets_[9].first + index;
    }

#if 0
    // auxiliary
    // [F^k is numeric] is equiv. to -F(k,t) for all t >= index_for_numerical
    void add_antecedent_for_numeric_feature(SAT::Implication &IP, int k) const {
        assert(0); // HOLA
        assert((0 <= k) && (k < K_));
        if( options_.log_feature_variables_ ) {
            for( int t = log_last_numerical_feature_; t < log_num_features_; ++t )
                IP.add_antecedent(-(1 + F(k, t)));
        } else {
            for( int j = S_.M().last_numerical_feature(); j < num_features_; ++j )
                IP.add_antecedent(-(1 + Fn(k, j)));
        }
    }
#endif

    // construction of base formulas
    void build_formulas_A() {
        // constraints for A formulas
        for( int j = 0; j < N_; ++j ) {
            for( int k = 0; k < K_; ++k ) {
                // if effect on F^k, it must be selected
                SAT::Implication *IP1 = new SAT::Implication;
                IP1->add_antecedent(1 + A(2, j, k));
                IP1->add_consequent(1 + A(1, j, k));
                add_implication(IP1);

                // if prec on F^k, it must be selected
                SAT::Implication *IP2 = new SAT::Implication;
                IP2->add_antecedent(1 + A(4, j, k));
                IP2->add_consequent(1 + A(3, j, k));
                add_implication(IP2);

                // if F^k is numeric and it is decremented, precondition must be value > 0;
                // that is, A(1,j,k) & -A(2,j,k) & [F^k is numeric] => A(4,j,k)
                SAT::Implication *IP3 = new SAT::Implication;
                IP3->add_antecedent(1 + A(1, j, k));
                IP3->add_antecedent(-(1 + A(2, j, k)));
                IP3->add_antecedent(1 + n(k));
                IP3->add_consequent(1 + A(4, j, k));
                add_implication(IP3);
            }
        }

        // some feature must be affected (selected)
        for( int j = 0; j < N_; ++j ) {
            SAT::Implication *IP = new SAT::Implication;
            for( int k = 0; k < K_; ++k )
                IP->add_consequent(1 + A(1, j, k));
            add_implication(IP);
        }

#if 0 // CHECK
        // selected action must be applicable at some state
        for( int j = 0; j < N_; ++j ) {
            SAT::Implication *IP = new SAT::Implication;
            for( int i = 0; i < num_states_; ++i )
                IP->add_consequent(1 + APP(i, j));
            add_implication(IP);
        }
#endif

        // ordering: if first affected feature by A^j is F^k, then A^{j-1} must affect some feature in {F^1, ..., F^k}
        if( options_.ordering_ && !options_.add_actions_ ) {
            for( int j = 1; j < N_; ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < k; ++t )
                        IP->add_antecedent(-(1 + A(1, j, t)));
                    IP->add_antecedent(1 + A(1, j, k));
                    for( int t = 0; t <= k; ++t )
                        IP->add_consequent(1 + A(1, j - 1, t));
                    add_implication(IP);
                }
            }
        }
    }

    void build_formulas_F() {
        if( options_.log_feature_variables_ ) {
            // constraints: dummy features are not selected
            for( int k = 0; k < K_; ++k ) {
                for( int index = S_.M().last_numerical_feature(); index < S_.M().first_boolean_feature(); ++index ) {
                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_antecedent(SAT::bit(t, index) ? 1 + F(k, t) : -(1 + F(k, t)));
                    add_implication(IP);
                }
            }

            // constraints: inexistent features are not selected
            for( int k = 0; k < K_; ++k ) {
                for( int index = num_features_; index < (1 << log_num_features_); ++index ) {
                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_antecedent(SAT::bit(t, index) ? 1 + F(k, t) : -(1 + F(k, t)));
                    add_implication(IP);
                }
            }
        } else {
            // F^k should be some feature
            for( int k = 0; k < K_; ++k ) {
                SAT::Implication *IP = new SAT::Implication;
                for( int j = 0; j < num_features_; ++j )
                    IP->add_consequent(1 + Fn(k, j));
                add_implication(IP);
            }

            // F^k cannot map to two different features
            for( int k = 0; k < K_; ++k ) {
                for( int i = 0; i < num_features_; ++i ) {
                    for( int j = i + 1; j < num_features_; ++j ) {
                        SAT::Implication *IP = new SAT::Implication;
                        IP->add_consequent(-(1 + Fn(k, i)));
                        IP->add_consequent(-(1 + Fn(k, j)));
                        add_implication(IP);
                    }
                }
            }

            // constraints: dummy features are not selected
            for( int k = 0; k < K_; ++k ) {
                for( int j = S_.M().last_numerical_feature(); j < S_.M().first_boolean_feature(); ++j ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + Fn(k, j));
                    add_implication(IP);
                }
            }
        }

        // ordering of selected features F^0 < F^0 < ... < F^{K-1}
        build_formulas_ordering();
    }

    void build_n(int k) {
        if( options_.log_feature_variables_ ) {
            SAT::Implication *IP1 = new SAT::Implication;
            for( int t = log_last_numerical_feature_; t < log_num_features_; ++t )
                IP1->add_antecedent(-(1 + F(k, t)));
            IP1->add_consequent(1 + n(k));
            add_implication(IP1);

            for( int t = log_last_numerical_feature_; t < log_num_features_; ++t ) {
                SAT::Implication *IP2 = new SAT::Implication;
                IP2->add_antecedent(1 + F(k, t));
                IP2->add_consequent(-(1 + n(k)));
                add_implication(IP2);
            }
        } else {
            for( int j = 0; j < num_features_; ++j ) {
                SAT::Implication *IP = new SAT::Implication;
                IP->add_antecedent(1 + Fn(k, j));
                IP->add_consequent(j < S_.M().last_numerical_feature() ? 1 + n(k) : -(1 + n(k)));
                add_implication(IP);
            }
        }
    }
    void build_formulas_n() {
        for( int k = 0; k < K_; ++k )
            build_n(k);
    }

    void build_phi(int k, int i, bool star) {
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));

        if( options_.log_feature_variables_ ) {
            // forward implications
            for( int l = 0; l < num_features_; ++l ) {
                if( (!star && (S_.M()(i, l) == 0)) || (star && (S_.M()(i, l) > 0)) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + phi(k, i, star));
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_consequent(SAT::bit(t, l) ? -(1 + F(k, t)) : (1 + F(k, t)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int l = 0; l < num_features_; ++l ) {
                if( (!star && (S_.M()(i, l) > 0)) || (star && (S_.M()(i, l) == 0)) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_consequent(1 + phi(k, i, star));
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_antecedent(SAT::bit(t, l) ? (1 + F(k, t)) : -(1 + F(k, t)));
                    add_implication(IP);
                }
            }
        } else {
            // forward implications
            for( int l = 0; l < num_features_; ++l ) {
                if( (!star && (S_.M()(i, l) == 0)) || (star && (S_.M()(i, l) > 0)) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + phi(k, i, star));
                    IP->add_consequent(-(1 + Fn(k, l)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int l = 0; l < num_features_; ++l ) {
                if( (!star && (S_.M()(i, l) > 0)) || (star && (S_.M()(i, l) == 0)) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_consequent(1 + phi(k, i, star));
                    IP->add_antecedent(1 + Fn(k, l));
                    add_implication(IP);
                }
            }
        }
    }
    void build_formulas_phi() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                build_phi(k, i, false);
                build_phi(k, i, true);
            }
        }
    }

    void build_B(int t, int j, int i, int k) { // t in {1,...,2}, j in {0, ..., N-1}, i in {0, ..., states-1}, k in {0, ..., K-1}
        assert((1 <= t) && (t < 3));
        assert((0 <= j) && (j < N_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= k) && (k < K_));
        // B(1,j,i,k) = A(4,j,k) -> phi(k,i)
        // B(2,j,i,k) = A(3,j,k) & -A(4,j,k) -> phi*(k,i)

        bool star = t == 1;

        // forward implication: one
        SAT::Implication *IP1 = new SAT::Implication;
        IP1->add_antecedent(1 + B(t, j, i, k));
        if( star ) IP1->add_antecedent(1 + A(3, j, k));
        IP1->add_antecedent(star ? -(1 + A(4, j, k)) : 1 + A(4, j, k));
        IP1->add_consequent(1 + phi(k, i, star));
        add_implication(IP1);

        // backward implications: two or three
        SAT::Implication *IP2 = new SAT::Implication;
        IP2->add_antecedent(star ? 1 + A(4, j, k) : -(1 + A(4, j, k)));
        IP2->add_consequent(1 + B(t, j, i, k));
        add_implication(IP2);

        SAT::Implication *IP3 = new SAT::Implication;
        IP3->add_antecedent(1 + phi(k, i, star));
        IP3->add_consequent(1 + B(t, j, i, k));
        add_implication(IP3);

        if( star ) {
            SAT::Implication *IP4 = new SAT::Implication;
            IP4->add_antecedent(-(1 + A(3, j, k)));
            IP4->add_consequent(1 + B(t, j, i, k));
            add_implication(IP4);
        }
    }
    void build_formulas_B() {
        for( int t = 1; t < 3; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int i = 0; i < num_states_; ++i ) {
                    for( int k = 0; k < K_; ++k )
                        build_B(t, j, i, k);
                }
            }
        }
    }

    void build_RES(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int state_l = S_.T().transition_state(i, l);

        // forward implications: 4 of them
        for( int k = 0; k < K_; ++k ) {
            // RES(i,l,j) & A(2,j,k) -> phi(k,l)
            SAT::Implication *IP1 = new SAT::Implication;
            IP1->add_antecedent(1 + RES(i, l, j));
            IP1->add_antecedent(1 + A(2, j, k));
            IP1->add_consequent(1 + phi(k, state_l, false));
            add_implication(IP1);

            // RES(i,l,j) & A(1,j,k) & -A(2,j,k) & -[F^k is numeric] -> phi*(k,l)
            SAT::Implication *IP2 = new SAT::Implication;
            IP2->add_antecedent(1 + RES(i, l, j));
            IP2->add_antecedent(1 + A(1, j, k));
            IP2->add_antecedent(-(1 + A(2, j, k)));
            IP2->add_antecedent(-(1 + n(k)));
            IP2->add_consequent(1 + phi(k, state_l, true));
            add_implication(IP2);

            // RES(i,l,j) & -A(1,j,k) & phi(k,i) -> phi(k,l)
            SAT::Implication *IP3 = new SAT::Implication;
            IP3->add_antecedent(1 + RES(i, l, j));
            IP3->add_antecedent(-(1 + A(1, j, k)));
            IP3->add_antecedent(1 + phi(k, i, false));
            IP3->add_consequent(1 + phi(k, state_l, false));
            add_implication(IP3);

            // RES(i,l,j) & -A(1,j,k) & phi*(k,i) -> phi*(k,l)
            SAT::Implication *IP4 = new SAT::Implication;
            IP4->add_antecedent(1 + RES(i, l, j));
            IP4->add_antecedent(-(1 + A(1, j, k)));
            IP4->add_antecedent(1 + phi(k, i, true));
            IP4->add_consequent(1 + phi(k, state_l, true));
            add_implication(IP4);
        }
    }
    void build_formulas_RES() {
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l )
                    build_RES(i, l, j);
            }
        }
    }

    void build_INC(int k, int i, int l) { // k in {0, ..., K-1}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        int state_l = S_.T().transition_state(i, l);

        if( options_.log_feature_variables_ ) {
            // forward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) >= S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + INC(k, i, l));
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_consequent(SAT::bit(t, j) ? -(1 + F(k, t)) : (1 + F(k, t)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) < S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_antecedent(SAT::bit(t, j) ? 1 + F(k, t) : -(1 + F(k, t)));
                    IP->add_consequent(1 + INC(k, i, l));
                    add_implication(IP);
                }
            }
        } else {
            // forward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) >= S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + INC(k, i, l));
                    IP->add_consequent(-(1 + Fn(k, j)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) < S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + Fn(k, j));
                    IP->add_consequent(1 + INC(k, i, l));
                    add_implication(IP);
                }
            }
        }
    }
    void build_formulas_INC() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l )
                    build_INC(k, i, l);
            }
        }
    }

    void build_DEC(int k, int i, int l) { // k in {0, ..., K-1}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        int state_l = S_.T().transition_state(i, l);

        if( options_.log_feature_variables_ ) {
            // forward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) <= S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + DEC(k, i, l));
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_consequent(SAT::bit(t, j) ? -(1 + F(k, t)) : (1 + F(k, t)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) > S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_antecedent(SAT::bit(t, j) ? 1 + F(k, t) : -(1 + F(k, t)));
                    IP->add_consequent(1 + DEC(k, i, l));
                    add_implication(IP);
                }
            }
        } else {
            // forward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) <= S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + DEC(k, i, l));
                    IP->add_consequent(-(1 + Fn(k, j)));
                    add_implication(IP);
                }
            }

            // backward implications
            for( int j = 0; j < num_features_; ++j ) {
                if( S_.M()(i, j) > S_.M()(state_l, j) ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + Fn(k, j));
                    IP->add_consequent(1 + DEC(k, i, l));
                    add_implication(IP);
                }
            }
        }
    }
    void build_formulas_DEC() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < S_.T().num_transitions(i); ++l )
                    build_DEC(k, i, l);
            }
        }
    }

    void build_MATCH(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        //int state_l = S_.T().transition_state(i, l); // unused

        // forward: implication for RES
        SAT::Implication *IP = new SAT::Implication;
        IP->add_antecedent(1 + MATCH(i, l, j));
        IP->add_consequent(1 + RES(i, l, j));
        add_implication(IP);

        // forward: increments for INCs
        for( int k = 0; k < K_; ++k ) {
            SAT::Implication *IP1 = new SAT::Implication;
            IP1->add_antecedent(1 + MATCH(i, l, j));
            IP1->add_antecedent(1 + A(2, j, k));
            IP1->add_antecedent(1 + n(k));
            IP1->add_consequent(1 + INC(k, i, l));
            add_implication(IP1);

            SAT::Implication *IP2 = new SAT::Implication;
            IP2->add_antecedent(1 + MATCH(i, l, j));
            IP2->add_antecedent(1 + INC(k, i, l));
            IP2->add_antecedent(1 + n(k));
            IP2->add_consequent(1 + A(2, j, k));
            add_implication(IP2);
        }

        // forward: increments for DECs
        for( int k = 0; k < K_; ++k ) {
            SAT::Implication *IP1 = new SAT::Implication;
            IP1->add_antecedent(1 + MATCH(i, l, j));
            IP1->add_antecedent(1 + A(1, j, k));
            IP1->add_antecedent(-(1 + A(2, j, k)));
            IP1->add_antecedent(1 + n(k));
            IP1->add_consequent(1 + DEC(k, i, l));
            add_implication(IP1);

            SAT::Implication *IP2 = new SAT::Implication;
            IP2->add_antecedent(1 + MATCH(i, l, j));
            IP2->add_antecedent(1 + DEC(k, i, l));
            IP2->add_antecedent(1 + n(k));
            IP2->add_consequent(1 + A(1, j, k));
            add_implication(IP2);

            SAT::Implication *IP3 = new SAT::Implication;
            IP3->add_antecedent(1 + MATCH(i, l, j));
            IP3->add_antecedent(1 + DEC(k, i, l));
            IP3->add_antecedent(1 + n(k));
            IP3->add_consequent(-(1 + A(2, j, k)));
            add_implication(IP3);
        }
    }
    void build_formulas_MATCH() {
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                for( int j = 0; j < N_; ++j )
                    build_MATCH(i, l, j);
            }
        }
    }

    void build_APPMATCH(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < S_.T().num_transitions(i)));
        assert((0 <= j) && (j < N_));
        //int state_l = S_.T().transition_state(i, l); // unused

        // forward implications
        SAT::Implication *IP = new SAT::Implication;
        IP->add_antecedent(1 + APPMATCH(i, l, j));
        IP->add_consequent(1 + MATCH(i, l, j));
        add_implication(IP);

        for( int k = 0; k < K_; ++k ) {
            for( int t = 1; t < 3; ++t ) {
                SAT::Implication *IP = new SAT::Implication;
                IP->add_antecedent(1 + APPMATCH(i, l, j));
                IP->add_consequent(1 + B(t, j, i, k));
                add_implication(IP);
            }
        }
    }
    void build_formulas_APPMATCH() {
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                for( int j = 0; j < N_; ++j )
                    build_APPMATCH(i, l, j);
            }
        }
    }

    // this formula says F^k < F^{k+1}
    void build_df(int k, int t) { // k in {0, ..., K-2}, t in {0, ..., log(F)}
        assert((0 <= k) && (k + 1 < K_));
        assert((0 <= t) && (t < log_num_features_));

        SAT::Implication *IP1 = new SAT::Implication;
        IP1->add_antecedent(1 + df(k, t));
        IP1->add_consequent(-(1 + F(k, t)));
        add_implication(IP1);

        SAT::Implication *IP2 = new SAT::Implication;
        IP2->add_antecedent(1 + df(k, t));
        IP2->add_consequent(1 + F(k + 1, t));
        add_implication(IP2);

        for( int j = t + 1; j < log_num_features_; ++j ) {
            SAT::Implication *IQ1 = new SAT::Implication;
            IQ1->add_antecedent(1 + df(k, t));
            IQ1->add_antecedent(1 + F(k, j));
            IQ1->add_consequent(1 + F(k + 1, j));
            add_implication(IQ1);

            SAT::Implication *IQ2 = new SAT::Implication;
            IQ2->add_antecedent(1 + df(k, t));
            IQ2->add_antecedent(-(1 + F(k, j)));
            IQ2->add_consequent(-(1 + F(k + 1, j)));
            add_implication(IQ2);
        }
    }
    void build_formulas_ordering() {
        if( do_ordering() ) {
            if( options_.log_feature_variables_ ) {
                for( int k = 0; k + 1 < K_; ++k ) {
                    for( int t = 0; t < log_num_features_; ++t )
                        build_df(k, t);

                    SAT::Implication *IP = new SAT::Implication;
                    for( int t = 0; t < log_num_features_; ++t )
                        IP->add_consequent(1 + df(k, t));
                    add_implication(IP);
                }
            } else {
                for( int k = 0; k + 1 < K_; ++k ) {
                    for( int i = 0; i < num_features_; ++i ) {
                        SAT::Implication *IP = new SAT::Implication;
                        IP->add_antecedent(1 + Fn(k, i));
                        for( int j = i + 1; j < num_features_; ++j )
                            IP->add_consequent(1 + Fn(k + 1, j));
                        add_implication(IP);
                    }
                }
            }
        }
    }

    virtual void build_base() {
        imp_offsets_.push_back(make_pair(0, "A"));
        add_comment("Theory for formulas A(t,j,k)");
        build_formulas_A();
        cout << "A: #implications=" << implications_.size() - imp_offsets_.back().first << " (3NK + N + (N-1) x K = N(1 + 4K) - K)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "F"));
        add_comment("Theory for formulas F(k,t)");
        build_formulas_F();
        cout << "F: #implications=" << implications_.size() - imp_offsets_.back().first << " (K x (D+U) + (K-1)) where D=#dummy-features and U=#unused-features" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
        //
        imp_offsets_.push_back(make_pair(implications_.size(), "n"));
        add_comment("Theory for formulas n(k)");
        build_formulas_n();
        cout << "n: #implications=" << implications_.size() - imp_offsets_.back().first << " (?)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "phi"));
        add_comment("Theory for formulas phi(k,i)");
        build_formulas_phi();
        cout << "phi: #implications=" << implications_.size() - imp_offsets_.back().first << " (2SK x F)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "B"));
        add_comment("Theory for formulas B(t,j,i,k)");
        build_formulas_B();
        cout << "B: #implications=" << implications_.size() - imp_offsets_.back().first << " (7SKN)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "RES"));
        add_comment("Theory for formulas RES(i,l,j)");
        build_formulas_RES();
        cout << "RES: #implications=" << implications_.size() - imp_offsets_.back().first << " (3TNK + TNK x Delta(L)) where Delta(L)=log2(#features)-log2(last_numerical_feature)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "INC"));
        add_comment("Theory for formulas INC(k,i,l)");
        build_formulas_INC();
        cout << "INC: #implications=" << implications_.size() - imp_offsets_.back().first << " (TKF)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "DEC"));
        add_comment("Theory for formulas DEC(k,i,l)");
        build_formulas_DEC();
        cout << "DEC: #implications=" << implications_.size() - imp_offsets_.back().first << " (TKF)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "MATCH"));
        add_comment("Theory for formulas MATCH(i,l,j)");
        build_formulas_MATCH();
        cout << "MATCH: #implications=" << implications_.size() - imp_offsets_.back().first << " (TN + 5TNK)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "APPMATCH"));
        add_comment("Theory for formulas APPMATCH(i,l,j)");
        build_formulas_APPMATCH();
        cout << "APPMATCH: #implications=" << implications_.size() - imp_offsets_.back().first << " (TN + 2TNK)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

#if 0
        // CHECK: for some reason, it is more efficient to emit this formula
        // CHECK: when emiting F formulas
        imp_offsets_.push_back(make_pair(implications_.size(), "ordering"));
        add_comment("Theory for ordering formulas");
        build_formulas_ordering();
        cout << "ordering: #implications=" << implications_.size() - imp_offsets_.back().first << " ((K-1) x (2L + L(L-1)))" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
#endif
    }

    // features added (forced)
    // j-th forced feature is feature whose index is feature_index
    void force_feature(int j, int feature_index) {
        if( do_ordering() ) {
            // some F^k is this feature
            SAT::Implication *IP = new SAT::Implication;
            for( int k = 0; k < K_; ++k )
                IP->add_consequent(1 + MappedBy(j, k));
            add_implication(IP);

            if( options_.log_feature_variables_ ) {
                // if F^k is feature, set bits F(k, t) accordingly
                for( int k = 0; k < K_; ++k ) {
                    for( int t = 0; t < log_num_features_; ++t ) {
                        SAT::Implication *IP = new SAT::Implication;
                        IP->add_antecedent(1 + MappedBy(j, k));
                        IP->add_consequent(SAT::bit(t, feature_index) ? (1 + F(k, t)) : -(1 + F(k, t)));
                        add_implication(IP);
                    }
                }
            } else {
                for( int k = 0; k < K_; ++k ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_antecedent(1 + MappedBy(j, k));
                    IP->add_consequent(1 + Fn(k, feature_index));
                    add_implication(IP);
                }
            }
        } else {
            if( options_.log_feature_variables_ ) {
                // F^j is j-th forced feature, set bits F(j, t) accordingly
                for( int t = 0; t < log_num_features_; ++t ) {
                    SAT::Implication *IP = new SAT::Implication;
                    IP->add_consequent(SAT::bit(t, feature_index) ? 1 + F(j, t) : -(1 + F(j, t)));
                    add_implication(IP);
                }
            } else {
                SAT::Implication *IP = new SAT::Implication;
                IP->add_consequent(1 + Fn(j, feature_index));
                add_implication(IP);
            }
        }
    }
    void force_feature(int j, const string &feature) {
        int feature_index = S_.M().feature(feature);
        if( feature_index == -1 )
            cout << "error: inexistent feature '" << feature << "'" << endl;
        else
            force_feature(j, feature_index);
    }
    void add_implications_for_features() {
        assert(forced_features_.size() <= K_);
        if( !forced_features_.empty() ) {
            imp_offsets_.push_back(make_pair(implications_.size(), "forced-features"));
            add_comment("Implications for added (forced) features");
        }

        int j = 0;
        for( Items<string>::const_iterator it = forced_features_.begin(); it != forced_features_.end(); ++it, ++j )
            force_feature(j, *it);

        if( do_ordering() )
            cout << "features: #implications=" << implications_.size() - imp_offsets_.back().first << " (M x (1 + K x log2(#features)))" << endl;
        else
            cout << "features: #implications=" << implications_.size() - imp_offsets_.back().first << " (M x log2(#features))" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }
    void add_variables_for_features() {
        if( do_ordering() ) {
            var_offsets_.push_back(make_pair(variables_.size(), "MappedBy"));
            for( int j = 0; j < int(forced_features_.size()); ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    int index = variables_.size();
                    string name = string("MappedBy(") + to_string(j) + "," + to_string(k) + ")";
                    variables_.push_back(new SAT::Var(index, name));
                    assert(index == MappedBy(j, k));
                }
            }
            cout << "MappedBy: #variables=" << variables_.size() - var_offsets_.back().first << " (MK), offset=" << var_offsets_.back().first << endl;
        }
    }
    void read_features(const string &prefix) {
        string features_filename = filename(prefix, K_, N_, "_features.dat", false);
        cout << "reading '" << features_filename << "' ... " << flush;
        ifstream ifs(features_filename.c_str());
        if( !ifs.fail() ) {
            const Items<string> *features = Items<string>::read_dump(ifs);
            ifs.close();
            forced_features_ = *features;
            delete features;
        } else {
            cout << "error: opening file '" << features_filename << "'" << endl;
        }
    }

    // CHECK: clauses to separate explicit goals from non-goals
    // Suggestion:
    //   1. Define variable not-sim(g,n,k) for goal state g and non-goal state n using feature k
    //   2. Require for each g, n, that there is f such that not-sim(g,n,k). Each requirement is
    //      a disjunction over k for not-sim(g,n,k)
    //   3. not-sim(g,n,k) need two implications:
    //      phi(g,f) & phi*(g,f) -> not-sim(g,n,k)
    //      phi*(g,f) & phi(g,f) -> not-sim(g,n,k)
    //   In total need G x N x L vars and G x N x (1 + 2K) implications where G=#goals and N=#non-goals
    void add_implications_for_goals() {
        assert(forced_goals_.empty());
    }
    void add_variables_for_goals() {
        //CHECK assert(0);
    }
    void read_goals(const string &prefix) {
        string goals_filename = filename(prefix, K_, N_, "_goals.dat", false);
        cout << "reading '" << goals_filename << "' ... " << flush;
        ifstream ifs(goals_filename.c_str());
        if( !ifs.fail() ) {
            const Items<int> *goals = Items<int>::read_dump(ifs);
            ifs.close();
            forced_goals_ = *goals;
            delete goals;
        } else {
            cout << "error: opening file '" << goals_filename << "'" << endl;
        }
    }

    // actions added (forced)
    void force_action(int j, const AbstractAction &action) {
        add_comment(string("Forced action ") + action.name_);
        for( int k = 0; k < K_; ++k ) {
            SAT::Implication *IP1 = new SAT::Implication;
            IP1->add_consequent(action.selected_effect_ & (1 << k) ? 1 + A(1, j, k) : -(1 + A(1, j, k)));
            add_implication(IP1);
            SAT::Implication *IP2 = new SAT::Implication;
            IP2->add_consequent(action.effect_ & (1 << k) ? 1 + A(2, j, k) : -(1 + A(2, j, k)));
            add_implication(IP2);
            SAT::Implication *IP3 = new SAT::Implication;
            IP3->add_consequent(action.selected_precondition_ & (1 << k) ? 1 + A(3, j, k) : -(1 + A(3, j, k)));
            add_implication(IP3);
            SAT::Implication *IP4 = new SAT::Implication;
            IP4->add_consequent(action.precondition_ & (1 << k) ? 1 + A(4, j, k) : -(1 + A(4, j, k)));
            add_implication(IP4);
        }
    }
    void add_implications_for_actions() {
        assert(forced_actions_.size() <= N_);

        imp_offsets_.push_back(make_pair(implications_.size(), "force-actions"));
        add_comment("Theory for forced actions");

        int j = 0;
        for( Items<AbstractAction>::const_iterator it = forced_actions_.begin(); it != forced_actions_.end(); ++it, ++j )
            force_action(j, *it);

        cout << "actions: #implications=" << implications_.size() - imp_offsets_.back().first << " (?)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }
    void read_actions(const string &prefix); // defined below

    void enforce_soundness() {
        imp_offsets_.push_back(make_pair(implications_.size(), "soundness"));
        add_comment("Theory for enforcing soundness");
        for( int i = 0; i < num_states_; ++i ) {
            for( int j = 0; j < N_; ++j ) {
                SAT::Implication *IP = new SAT::Implication;

                // antecedent
                for( int k = 0; k < K_; ++k ) {
                    for( int t = 1; t < 3; ++t )
                        IP->add_antecedent(1 + B(t, j, i, k));
                }

                // consequent
                for( int l = 0; l < S_.T().num_transitions(i); ++l )
                    IP->add_consequent(1 + MATCH(i, l, j));

                // add implications
                add_implication(IP);
            }
        }
        cout << "soundness: #implications=" << implications_.size() - imp_offsets_.back().first << " (SN)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }

    void enforce_completeness() {
        imp_offsets_.push_back(make_pair(implications_.size(), "completeness"));
        add_comment("Theory for enforcing completeness");
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < S_.T().num_transitions(i); ++l ) {
                SAT::Implication *IP = new SAT::Implication;
                for( int j = 0; j < N_; ++j )
                    IP->add_consequent(1 + APPMATCH(i, l, j));
                add_implication(IP);
            }
        }
        cout << "completeness: #implications=" << implications_.size() - imp_offsets_.back().first << " (T)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }

  public:
    // print abstraction coded in model
    virtual void decode_model(std::ostream &os) const {
        decode_model("<empty>", os, true);
    }
    void decode_model(const string &name, ostream &os, bool qnp_format) const {
        assert(satisfiable_ && (model_.size() == variables_.size()));

        // name
        os << name << endl;
        if( qnp_format )
            os << K_;
        else
            os << "Features:" << endl;

        // features
        for( int j = 0; j < K_; ++j ) {
            int feature_index = index_for_selected_feature(model_, j);
            assert(feature_index < num_features_);
            if( qnp_format )
                os << " " << S_.M().feature(feature_index) << " " << S_.M().numerical(feature_index);
            else
                os << "F" << j << ": " << S_.M().feature(feature_index) << endl;
        }
        if( qnp_format ) os << endl;

        // abstract actions
        set<AbstractAction> unique;
        for( int j = 0; j < N_; ++j ) {
            AbstractAction act(string("decoded-") + to_string(j));
            for( int k = 0; k < K_; ++k ) {
                act.selected_precondition_ += model_[A(3, j, k)] << k;
                act.precondition_ += model_[A(4, j, k)] << k;
                act.selected_effect_ += model_[A(1, j, k)] << k;
                act.effect_ += model_[A(2, j, k)] << k;
            }
            unique.insert(act);
        }

        if( qnp_format ) {
            os << unique.size() << endl;
            for( set<AbstractAction>::const_iterator it = unique.begin(); it != unique.end(); ++it )
                it->print(os, *this, model_, qnp_format);
        } else {
            os << "Abstract actions (#total=" << N_ << ", #unique=" << unique.size() << "):" << endl;
            for( set<AbstractAction>::const_iterator it = unique.begin(); it != unique.end(); ++it )
                it->print(os, *this, model_, qnp_format);
        }
    }
};

// static member
map<string, int> *Theory::AbstractAction::feature_map_ = nullptr;

istream& operator>>(istream &is, Theory::AbstractAction &act) {
    string name;
    is >> name;
    act.name_ = name;

    int nprec, neffect;
    for( is >> nprec; nprec > 0; --nprec ) {
        string feature;
        int value;
        is >> feature >> value;
        act.add_precondition(feature, value);
    }
    for( is >> neffect; neffect > 0; --neffect ) {
        string feature;
        int value;
        is >> feature >> value;
        act.add_effect(feature, value);
    }
  
    return is;
}

void Theory::read_actions(const string &prefix) {
    // construct and set feature map
    int j = 0;
    map<string, int> feature_map;
    for( Items<string>::const_iterator it = forced_features_.begin(); it != forced_features_.end(); ++it, ++j )
        feature_map.insert(make_pair(*it, j));
    Theory::AbstractAction::feature_map_ = &feature_map;

    // read file
    string actions_filename = filename(prefix, K_, N_, "_actions.dat", false);
    cout << "reading '" << actions_filename << "' ... " << flush;
    ifstream ifs(actions_filename.c_str());
    if( !ifs.fail() ) {
        const Items<AbstractAction> *actions = Items<AbstractAction>::read_dump(ifs);
        ifs.close();
        forced_actions_ = *actions;
        delete actions;
    } else {
        cout << "error: opening file '" << actions_filename << "'" << endl;
    }

    // reset feature map
    Theory::AbstractAction::feature_map_ = nullptr;
}

void usage(ostream &os, const string &name) {
    cout << endl
         << "usage: " << name << " [--add-actions] [--add-features] [--decode] [--explicit-goals] [--naive-feature-variables] [--only-soundness] [--ordering-clauses] <prefix> <K> <N>" << endl
         << endl
         << "where" << endl
         << "    <prefix> is prefix for all files" << endl
         << "    K is numer of features to select" << endl
         << "    N is number of abstract actions." << endl
         << endl
         << "For the options," << endl
         << "    --add-actions to force inclusion of actions in '<prefix>_actions.dat'" << endl
         << "    --add-features to force inclusion of features in '<prefix>_features.dat'" << endl
         << "    --decode to decode model in '<prefix>_<K>_<N>_model.cnf' found by minisat" << endl
         << "    --explicit-goals to separate goals from non-goals using goals in '<prefix>_goals.dat'" << endl
         << "    --naive-feature-variables to generate naive encoding of feature variables rather than log encoding" << endl
         << "    --only-soundness to generate clauses that only enforce soundness instead of soundness+completeness" << endl
         << "    --ordering-clauses to generate clauses to enforce ordering of selected features and abstract actions" << endl
         << endl
         << "A feature file consists of number <n> of features in file, followed by space" << endl
         << "followed by space-separated list of features names (the names need to match those" << endl
         << "in matrix file." << endl
         << endl
         ;
}

int main(int argc, const char **argv) {
    // print call
    cout << "call:";
    for( int i = 0; i < argc; ++i )
        cout << " " << argv[i];
    cout << endl;

    // read executable name
    string name(*argv);

    // read options
    Theory::Options options;
    for( ++argv, --argc; (argc > 0) && (**argv == '-'); ++argv, --argc ) {
        if( string(*argv) == "--decode" ) {
            options.decode_ = true;
        } else if( string(*argv) == "--explicit-goals" ) {
            options.explicit_goals_ = true;
        } else if( string(*argv) == "--add-actions" ) {
            options.add_actions_ = true;
        } else if( string(*argv) == "--add-features" ) {
            options.add_features_ = true;
        } else if( string(*argv) == "--naive-feature-variables" ) {
            options.log_feature_variables_ = false;
        } else if( string(*argv) == "--only-soundness" ) {
            options.only_soundness_ = true;
        } else if( string(*argv) == "--ordering-clauses" ) {
            options.ordering_ = true;
        } else {
            cout << "error: unrecognized option '" << *argv << "'" << endl;
            usage(cout, name);
            exit(0);
        }
    }

    // check we have enough arguments
    if( argc < 3 ) {
        usage(cout, name);
        exit(0);
    }

    // read arguments
    options.prefix_ = argv[0];
    int K = atoi(argv[1]);
    int N = atoi(argv[2]);
    cout << "arguments: prefix=" << options.prefix_ << ", K=" << K << ", N=" << N << endl;

    // input filenames
    string matrix_filename = filename(options.prefix_, K, N, "_matrix.dat", false);
    string transitions_filename = filename(options.prefix_, K, N, "_transitions.dat", false);

    // read data
    Sample::Sample data;
    data.read(matrix_filename, transitions_filename);

    // create theory
    Theory Th(data, K, N, options);
    cout << "#variables=" << Th.num_variables() << endl;
    cout << "#implications=" << Th.num_implications() << endl;

    // decode
    if( options.decode_ ) {
        string model_filename = filename(options.prefix_, K, N, "_model.cnf");
        cout << "reading file '" << model_filename << "' ..." << flush;
        ifstream ifs(model_filename.c_str());
        if( !ifs.fail() ) {
            Th.read_minisat_output(ifs);
            ifs.close();
            cout << " done!" << endl;

            // output model in qnp format
            string qnp_filename = filename(options.prefix_, K, N, ".qnp");
            cout << "writing file '" << qnp_filename << "' ..." << flush;
            ofstream os(qnp_filename.c_str());
            Th.decode_model(model_filename, os, true);
            os.close();
            cout << " done!" << endl;
        } else {
            cout << "error: opening file '" << model_filename << "'" << endl;
        }
    } else {
        // output theory in human readable format
        //Th.print(cout);

        // output theory in SATLIB format
        string theory_filename = filename(options.prefix_, K, N, "_theory.cnf");
        cout << "writing file '" << theory_filename << "' ..." << flush;
        ofstream os(theory_filename.c_str());
        Th.dump(os);
        os.close();
        cout << " done!" << endl;
    }

    return 0;
}

