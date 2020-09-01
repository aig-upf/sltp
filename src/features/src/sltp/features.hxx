
#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>
#include <boost/bimap.hpp>

#include <sltp/base.hxx>
#include <sltp/utils.hxx>
#include <sltp/algorithms.hxx>

#include <ctime>

namespace Sample {
    class TransitionSample;
}


namespace SLTP::DL {

const unsigned PRIMITIVE_COMPLEXITY = 1;

using object_id_t = unsigned;
using predicate_id_t = unsigned;
using atom_id_t = unsigned;
using state_id_t = unsigned;

class Atom;
class Sample;
class Concept;
class Role;
class Feature;
class State;

using feature_cache_t = std::unordered_map<feature_sample_denotation_t, const Feature*, utils::container_hash<feature_sample_denotation_t> >;


//! Command-line option processing
struct Options {
    std::string workspace;
    int timeout;
    unsigned complexity_bound;
    unsigned dist_complexity_bound;
    unsigned cond_complexity_bound;
    bool comparison_features;
    bool generate_goal_concepts;
    bool print_denotations;
};

// We cache sample and state denotations. The latter are cached
// into a simple hash (i.e. unordered_set). The former are cached
// using two hash maps (i.e. unordered_map): one that maps sample
// denotations to concept names, and the other that maps concept
// names to sample denotations.
//
// We also cache atoms that are used to represent states in the sample.
class Cache {
public:
    struct cache_support_t {
        // cache for sample denotations
        bool operator()(const sample_denotation_t *d1, const sample_denotation_t *d2) const {
            assert((d1 != nullptr) && (d2 != nullptr));
            assert(d1->size() == d2->size()); // number of states is fixed in sample
            return *d1 == *d2;
        }

        size_t operator()(const sample_denotation_t *obj) const {
            assert(obj != nullptr);
            size_t hash = (*this)((*obj)[0]);
            for (auto elem:*obj)
                hash = hash ^ (*this)(elem);
            return hash;
        }

        // cache for state denotations
        bool operator()(const state_denotation_t *sd1, const state_denotation_t *sd2) const {
            assert((sd1 != nullptr) && (sd2 != nullptr));
            // dimension of sd1 and sd2 may not be equal since they may
            // be for states with different number of objects, or for
            // denotation of concept/roles
            return *sd1 == *sd2;
        }

        size_t operator()(const state_denotation_t *obj) const {
            assert(obj != nullptr);
            std::hash<std::vector<bool> > hasher;
            return hasher(*static_cast<const std::vector<bool>*>(obj));
        }
    };

    using cache1_t = std::unordered_map<const sample_denotation_t*, std::string, cache_support_t, cache_support_t>;
    using cache2_t = std::unordered_map<std::string, const sample_denotation_t*>;
    using cache3_t = std::unordered_set<const state_denotation_t*, cache_support_t, cache_support_t>;

protected:
    cache1_t cache1_;
    cache2_t cache2_;
    cache3_t cache3_;

public:
    Cache() = default;
    ~Cache() = default;

    // cache1: (full) sample denotations for concepts
    //
    // We maintain a hash of (full) sample denotations so that one copy
    // exists for each such denotations. This hash is implemented with
    // two tables, one that provides mappings from sample denotations
    // to concept names, and the other is the inverse
    //
    // sample denotation -> concept name
    const cache1_t& cache1() const {
        return cache1_;
    }

    const sample_denotation_t* find_sample_denotation(const sample_denotation_t &d) const {
        auto it = cache1_.find(&d);
        return it == cache1_.end() ? nullptr : it->first;
    }

    const sample_denotation_t* find_or_insert_sample_denotation(const sample_denotation_t &d, const std::string &name) {
        auto it = cache1_.find(&d);
        if( it == cache1_.end() ) {
            assert(cache2_.find(name) == cache2_.end());
            const sample_denotation_t *nd = new sample_denotation_t(d);
            cache1_.emplace(nd, name);
            cache2_.emplace(name, nd);
            return nd;
        } else {
            return it->first;
        }
    }
    const sample_denotation_t* find_or_insert_sample_denotation_by_name(const std::string &name, const sample_denotation_t &d) {
        auto it = cache2_.find(name);
        if( it == cache2_.end() ) {
            const sample_denotation_t *nd = new sample_denotation_t(d);
            cache2_.emplace(name, nd);
            return nd;
        } else {
            return it->second;
        }
    }

    // cache2: (full) sample denotations for concepts
    //
    // concept name -> sample denotation
    const sample_denotation_t* find_sample_denotation(const std::string &name) const {
        auto it = cache2_.find(name);
        return it == cache2_.end() ? nullptr : it->second;
    }

    // cache3: state denotations (bit vectors)
    //
    // We maintain a hash of these denotations such as to keep just one
    // copy of each such bitmap. A sample denotation is a vector of
    // state denotations. The same bitmap may appear associated with
    // different states in the same sample, or across sample denotations
    // for different concepts
    //
    // state denotation -> pointer
    const state_denotation_t* find_or_insert_state_denotation(const state_denotation_t &sd) {
        auto it = cache3_.find(&sd);
        if( it == cache3_.end() ) {
            const state_denotation_t *nsd = new state_denotation_t(sd);
            cache3_.insert(nsd);
            return nsd;
        } else {
            return *it;
        }
    }

    const state_denotation_t& retrieve_concept_denotation(const Concept &element, const State &state) const;
    const state_denotation_t& retrieve_role_denotation(const Role &element, const State &state) const;
};

// We represent states as subets of atoms
// An atom is a predicate and a vector of objects to ground the predicates.
class Object {
protected:
    const object_id_t id_;
    const std::string name_;

public:
    Object(unsigned id, std::string name) : id_(id), name_(std::move(name)) { }
    [[nodiscard]] int id() const {
        return id_;
    }
    [[nodiscard]] const std::string& as_str() const {
        return name_;
    }
    friend std::ostream& operator<<(std::ostream &os, const Object &obj) {
        return os << obj.as_str() << std::flush;
    }
};

struct Predicate {
    const predicate_id_t id_;
    const std::string name_;
    const int arity_;
    Predicate(unsigned id, std::string name, int arity)
            : id_(id), name_(std::move(name)), arity_(arity) {
    }
    [[nodiscard]] predicate_id_t id() const {
        return id_;
    }
    [[nodiscard]] const std::string& name() const {
        return name_;
    }
    [[nodiscard]] int arity() const {
        return arity_;
    }
    std::string as_str(const std::vector<object_id_t> *objects) const {
        std::string str = name_ + "(";
        if( objects == nullptr ) {
            for( int i = 0; i < arity_; ++i ) {
                str += std::string("x") + std::to_string(1 + i);
                if( 1 + i < arity_ ) str += ",";
            }
        } else {
            assert(objects->size() == arity_);
            for( int i = 0; i < arity_; ++i ) {
                //str += (*objects)[i]->as_str();
                str += std::to_string((*objects)[i]);
                if( 1 + i < arity_ ) str += ",";
            }
        }
        return str + ")";
    }
    friend std::ostream& operator<<(std::ostream &os, const Predicate &pred) {
        return os << pred.as_str(nullptr) << std::flush;
    }
};

class Atom {
protected:
    const predicate_id_t predicate_;
    const std::vector<object_id_t> objects_;

public:
    Atom(const predicate_id_t &predicate, std::vector<object_id_t> &&objects)
            : predicate_(predicate), objects_(std::move(objects)) {
    }

    [[nodiscard]] predicate_id_t pred_id() const {
        return predicate_;
    }
    [[nodiscard]] const std::vector<object_id_t>& objects() const {
        return objects_;
    }

    // Return the i-th object of the current atom
    [[nodiscard]] object_id_t object(int i) const {
        return objects_.at(i);
    }

    [[nodiscard]] bool is_instance(const Predicate &predicate) const {
        return predicate_ == predicate.id();
    }

    [[nodiscard]] std::vector<unsigned> data() const {
        std::vector<unsigned> res(1, predicate_);
        res.insert(res.end(), objects_.begin(), objects_.end());
        return res;
    }

    [[nodiscard]] std::string as_str(const Sample &sample) const;
};

// An instance stores information shared by the states that
// belong to the instance: objects and atoms mostly
class Instance {
public:
    // map from object name to object id in instance
//    using ObjectIndex = std::unordered_map<std::string, object_id_t>;
    using ObjectIndex = boost::bimap<std::string, object_id_t>;
    // map from atom of the form <pred_id, oid_1, ..., oid_n> to atom id in instance
    using AtomIndex = std::unordered_map<std::vector<unsigned>, atom_id_t, utils::container_hash<std::vector<unsigned> > >;

    const unsigned id;

protected:
    const std::vector<Object> objects_;
    const std::vector<Atom> atoms_;

    // mapping from object names to their ID in the sample
    ObjectIndex object_index_;

    // mapping from <predicate name, obj_name, ..., obj_name> to the ID of the corresponding GroundPredicate
    AtomIndex atom_index_;

public:
    Instance(unsigned id,
             std::vector<Object> &&objects,
             std::vector<Atom> &&atoms,
             ObjectIndex &&object_index,
             AtomIndex &&atom_index)
            :
        id(id),
        objects_(std::move(objects)),
        atoms_(std::move(atoms)),
        object_index_(std::move(object_index)),
        atom_index_(std::move(atom_index))
    {
    }

    Instance(const Instance& ins) = default;
    Instance(Instance &&ins) = default;
    ~Instance() = default;

    unsigned num_objects() const {
        return (unsigned) objects_.size();
    }
    int num_atoms() const {
        return atoms_.size();
    }
    const Atom& atom(unsigned id_) const {
        return atoms_.at(id_);
    }

    const std::vector<Atom>& atoms() const {
        return atoms_;
    }
    const ObjectIndex& object_index() const {
        return object_index_;
    }
    const AtomIndex& atom_index() const {
        return atom_index_;
    }
};

// A state is a collections of atoms
class State {
protected:
    const Instance &instance_;
    const state_id_t id_;
    std::vector<atom_id_t> atoms_;

public:
    explicit State(const Instance &instance, unsigned id, std::vector<atom_id_t> &&atoms)
            : instance_(instance), id_(id), atoms_(std::move(atoms)) {
    }
    State(const State &state) = default;
    State(State &&state) = default;

    [[nodiscard]] unsigned id() const {
        return id_;
    }
    [[nodiscard]] const std::vector<atom_id_t>& atoms() const {
        return atoms_;
    }
    [[nodiscard]] unsigned num_objects() const {
        return instance_.num_objects();
    }

    [[nodiscard]] const Instance& instance() const {
        return instance_;
    }

    [[nodiscard]] const Atom& atom(atom_id_t id) const {
        return instance_.atom(id);
    }
};

// A sample is a bunch of states and transitions among them. The
// sample contains the predicates used in the states, the objects,
// and the atoms
class Sample {
public:
    using PredicateIndex = std::unordered_map<std::string, predicate_id_t>;

protected:
    const std::vector<Predicate> predicates_;
    const std::vector<Instance> instances_;
    const std::vector<State> states_;

    // The IDs of predicates that are mentioned in the goal
    const std::vector<predicate_id_t> goal_predicates_;

    // mapping from predicate names to their ID in the sample
    PredicateIndex predicate_index_;

    Sample(std::vector<Predicate> &&predicates,
           std::vector<Instance> &&instances,
           std::vector<State> &&states,
           std::vector<predicate_id_t> &&goal_predicates,
           PredicateIndex &&predicate_index)
            : predicates_(std::move(predicates)),
              instances_(std::move(instances)),
              states_(std::move(states)),
              goal_predicates_(std::move(goal_predicates)),
              predicate_index_(std::move(predicate_index)) {
        std::cout << "SAMPLE:"
                  << " #predicates=" << predicates_.size()
                  << ", #instances=" << instances_.size()
                  << ", #states=" << states_.size()
                  << std::endl;
    }

public:
    Sample(const Sample &sample) = default;
    Sample(Sample &&sample) = default;
    ~Sample() = default;

    std::size_t num_predicates() const {
        return predicates_.size();
    }
    std::size_t num_states() const {
        return states_.size();
    }

    const std::vector<Predicate>& predicates() const {
        return predicates_;
    }
    const std::vector<predicate_id_t>& goal_predicates() const {
        return goal_predicates_;
    }
    const std::vector<State>& states() const {
        return states_;
    }

    const Predicate& predicate(predicate_id_t id) const {
        return predicates_.at(id);
    }
    const State& state(unsigned id) const {
        return states_.at(id);
    }
    const PredicateIndex& predicate_index() const {
        return predicate_index_;
    }

    // factory method - reads sample from serialized data
    static Sample read(std::istream &is);
};

inline std::string Atom::as_str(const Sample &sample) const {
    return sample.predicate(predicate_).as_str(&objects_);
}

//struct Denotation {
//    enum class denotation_t : bool { concept_denotation_t, role_denotation_t };
//    denotation_t type_;
//    std::vector<bool> values_;
//    Denotation(denotation_t type, size_t dimension)
//      : type_(type),
//        values_(std::vector<bool>(dimension, false)) {
//    }
//};
//
//// Denotation matrix is just a matrix of DL Denotations, each state being a row, each concept / role a column
//class DenotationMatrix {
//  protected:
//    std::vector<std::vector<Denotation> > data_;
//};
//
//class Model {
//};

class Base {
protected:
    int complexity_;

public:
    explicit Base(int complexity) : complexity_(complexity) { }
    [[nodiscard]] int complexity() const { return complexity_; }

    //virtual void denotation(Denotation &d) const = 0;
    virtual const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const = 0;
    virtual const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const = 0;
    [[nodiscard]] virtual std::string as_str() const = 0;

    [[nodiscard]] std::string as_str_with_complexity() const {
        return std::to_string(complexity_) + "." + as_str();
    }

    void force_complexity(int c) {
        complexity_ = c;
    }

    friend std::ostream& operator<<(std::ostream &os, const Base &base) {
        return os << base.as_str_with_complexity() << std::flush;
    }
};

class Concept : public Base {
public:
    explicit Concept(int complexity) : Base(complexity) { }
    virtual ~Concept() = default;
    [[nodiscard]] virtual const Concept* clone() const = 0;
};

class Role : public Base {
public:
    explicit Role(int complexity) : Base(complexity) { }
    virtual ~Role() = default;
    [[nodiscard]] virtual const Role* clone() const = 0;
};

class PrimitiveConcept : public Concept {
protected:
    const Predicate *predicate_;

public:
    explicit PrimitiveConcept(const Predicate *predicate) : Concept(PRIMITIVE_COMPLEXITY), predicate_(predicate) { }
    ~PrimitiveConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new PrimitiveConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            sample_denotation_t d;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                d.emplace_back(denotation(cache, sample, s));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(d, as_str()) : new sample_denotation_t(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        state_denotation_t sd(state.num_objects(), false);
        for( int i = 0; i < int(state.atoms().size()); ++i ) {
            atom_id_t id = state.atoms()[i];
            const Atom &atom = state.atom(id);
            if( atom.is_instance(*predicate_) ) {
                assert(atom.objects().size() == 1);
                object_id_t index = atom.object(0);
                assert(index < state.num_objects());
                sd[index] = true;
            }
        }
        return cache.find_or_insert_state_denotation(sd);
    }
    [[nodiscard]] std::string as_str() const override {
        return predicate_->name_;
    }
};


class NominalConcept : public Concept {
protected:
    const std::string& name_;

public:
    explicit NominalConcept(const std::string &name) : Concept(1), name_(name) {}
    ~NominalConcept() override = default;

    [[nodiscard]] const Concept* clone() const override {
        return new NominalConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            sample_denotation_t d;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                d.emplace_back(denotation(cache, sample, s));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(d, as_str()) : new sample_denotation_t(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        const Instance::ObjectIndex& oidx = state.instance().object_index();
        object_id_t id = oidx.left.at(name_);
        assert(id < state.num_objects());

        state_denotation_t sd(state.num_objects(), false);
        sd[id] = true;

        return cache.find_or_insert_state_denotation(sd);
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Nominal(") + name_ + ")";
    }
};

class UniversalConcept : public Concept {
public:
    UniversalConcept() : Concept(0) { }
    ~UniversalConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new UniversalConcept;
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            sample_denotation_t nd;
            nd.reserve(sample.num_states());
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                state_denotation_t sd(state.num_objects(), true);
                const state_denotation_t *cached_sd = cache.find_or_insert_state_denotation(sd);
                nd.emplace_back(cached_sd);
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: UniversalConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }
    [[nodiscard]] std::string as_str() const override {
        return "<universe>";
    }
};


class EmptyConcept : public Concept {
public:
    EmptyConcept() : Concept(0) { }
    ~EmptyConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new EmptyConcept;
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            sample_denotation_t nd;
            nd.reserve(sample.num_states());
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                state_denotation_t sd(state.num_objects(), false);
                const state_denotation_t *cached_sd = cache.find_or_insert_state_denotation(sd);
                nd.emplace_back(cached_sd);
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: EmptyConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }
    [[nodiscard]] std::string as_str() const override {
        return "<empty>";
    }
};

class AndConcept : public Concept {
protected:
    const Concept *concept1_;
    const Concept *concept2_;

public:
    AndConcept(const Concept *concept1, const Concept *concept2) :
            Concept(1 + concept1->complexity() + concept2->complexity()),
//        Concept(1 + concept1->complexity() * concept2->complexity()),
            concept1_(concept1),
            concept2_(concept2) {
    }
    ~AndConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new AndConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *d1 = cache.find_sample_denotation(concept1_->as_str());
            assert((d1 != nullptr) && (d1->size() == sample.num_states()));
            const sample_denotation_t *d2 = cache.find_sample_denotation(concept2_->as_str());
            assert((d2 != nullptr) && (d2->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sd1 = (*d1)[i];
                assert((sd1 != nullptr) && (sd1->size() == state.num_objects()));
                const state_denotation_t *sd2 = (*d2)[i];
                assert((sd2 != nullptr) && (sd2->size() == state.num_objects()));

                state_denotation_t nsd(state.num_objects(), false);
                for( int j = 0; j < state.num_objects(); ++j )
                    nsd[j] = (*sd1)[j] && (*sd2)[j];
                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: AndConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("And(") +  concept1_->as_str() + "," + concept2_->as_str() + ")";
    }
};

class NotConcept : public Concept {
protected:
    const Concept *concept_;

public:
    explicit NotConcept(const Concept *concept)
        : Concept(1 + concept->complexity()),
        concept_(concept)
    {}
    ~NotConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new NotConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sd = (*d)[i];
                assert((sd != nullptr) && (sd->size() == state.num_objects()));

                state_denotation_t nsd(state.num_objects(), false);
                for( int j = 0; j < state.num_objects(); ++j )
                    nsd[j] = !(*sd)[j];
                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: NotConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Not(") + concept_->as_str() + ")";
    }
};

class ExistsConcept : public Concept {
protected:
    const Concept *concept_;
    const Role *role_;

public:
    ExistsConcept(const Concept *concept, const Role *role)
        : Concept(1 + concept->complexity() + role->complexity()),
    concept_(concept),
    role_(role) {
    }
    ~ExistsConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new ExistsConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                unsigned m = state.num_objects();
                assert(state.id() == i);
                const state_denotation_t *c_den = (*d)[i];
                assert((c_den != nullptr) && (c_den->size() == m));
                const state_denotation_t *r_den = (*r)[i];
                assert((r_den != nullptr) && (r_den->size() == m*m));

                state_denotation_t nsd(m, false);
                for(unsigned x = 0; x < m; ++x) {

                    // x makes it into the denotation if there is an y such that y in c_den and (x,y) in r_den
                    for(unsigned y = 0; y < m; ++y) {
                        if((*c_den)[y]) {
                            auto x_y = x * m + y;
                            if ((*r_den)[x_y]) {
                                nsd[x] = true;
                                break;
                            }
                        }
                    }
                }

                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: ExistsConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Exists(") + role_->as_str() + "," + concept_->as_str() + ")";
    }
};

class ForallConcept : public Concept {
protected:
    const Concept *concept_;
    const Role *role_;

public:
    ForallConcept(const Concept *concept, const Role *role)
        : Concept(1 + concept->complexity() + role->complexity()),
    concept_(concept),
    role_(role) {
    }
    ~ForallConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new ForallConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                unsigned m = state.num_objects();
                assert(state.id() == i);
                const state_denotation_t *c_den = (*d)[i];
                assert((c_den != nullptr) && (c_den->size() == m));
                const state_denotation_t *r_den = (*r)[i];
                assert((r_den != nullptr) && (r_den->size() == m*m));

                state_denotation_t nsd(m, true);
                for(unsigned x = 0; x < m; ++x) {

                    // x does *not* make it into the denotation if there is an y
                    // such that y not in c_den and (x,y) in r_den
                    for(unsigned y = 0; y < m; ++y) {
                        if(!(*c_den)[y]) {
                            auto x_y = x * m + y;
                            if ((*r_den)[x_y]) {
                                nsd[x] = false;
                                break;
                            }
                        }
                    }
                }

                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: ForallConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Forall(") + role_->as_str() + "," + concept_->as_str() + ")";
    }
};


class EqualConcept : public Concept {
protected:
    const Role *r1_;
    const Role *r2_;

public:
    EqualConcept(const Role *r1, const Role *r2)
            : Concept(1 + r1->complexity() + r2->complexity()),
              r1_(r1),
              r2_(r2) {
    }
    ~EqualConcept() override = default;
    [[nodiscard]] const Concept* clone() const override { return new EqualConcept(*this); }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *r1_den = cache.find_sample_denotation(r1_->as_str());
            assert(r1_den && (r1_den->size() == sample.num_states()));
            const sample_denotation_t *r2_den = cache.find_sample_denotation(r2_->as_str());
            assert(r2_den && (r2_den->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                unsigned m = state.num_objects();
                assert(state.id() == i);
                const state_denotation_t *sd1 = (*r1_den)[i];
                assert((sd1 != nullptr) && (sd1->size() == m*m));
                const state_denotation_t *sd2 = (*r2_den)[i];
                assert((sd2 != nullptr) && (sd2->size() == m*m));

                state_denotation_t nsd(m, false);

                for(int x = 0; x < state.num_objects(); ++x) {
                    // If the set of y such that (x, y) in sd1 is equal to the set of z such that (x, z) in sd2,
                    // then x makes it into the denotation of this concept
                    bool in_denotation = true;
                    for (int z = 0; z < m; ++z) {
                        auto idx = x*m + z;
                        if ((*sd1)[idx] != (*sd2)[idx]) {
                            in_denotation = false;
                            break;
                        }
                    }

                    if (in_denotation) {
                        nsd[x] = true;
                    }
                }

                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }

    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: EqualConcept::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Equal(") + r1_->as_str() + "," + r2_->as_str() + ")";
    }
};


class PrimitiveRole : public Role {
protected:
    const Predicate *predicate_;

public:
    explicit PrimitiveRole(const Predicate *predicate) : Role(PRIMITIVE_COMPLEXITY), predicate_(predicate) { }
    ~PrimitiveRole() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new PrimitiveRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            sample_denotation_t nr;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                nr.emplace_back(denotation(cache, sample, s));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nr, as_str()) : new sample_denotation_t(nr);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        state_denotation_t sr(state.num_objects() * state.num_objects(), false);
        for( int i = 0; i < int(state.atoms().size()); ++i ) {
            atom_id_t id = state.atoms()[i];
            const Atom &atom = state.atom(id);
            if( atom.is_instance(*predicate_) ) {
                assert(atom.objects().size() == 2);
                int index = atom.object(0) * state.num_objects() + atom.object(1);
                assert(index < state.num_objects() * state.num_objects());
                sr[index] = true;
            }
        }
        return cache.find_or_insert_state_denotation(sr);
    }

    [[nodiscard]] const Predicate* predicate() const { return predicate_; }

    [[nodiscard]] std::string as_str() const override {
        return predicate_->name_;
    }
};

class PlusRole : public Role {
protected:
    const Role *role_;

public:
    explicit PlusRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    ~PlusRole() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new PlusRole(*this);
    }

    // apply Johnson's algorithm for transitive closure
    static void transitive_closure(int num_objects, state_denotation_t &sd) {
        // create adjacency lists
        std::vector<std::vector<int> > adj(num_objects);
        for( int i = 0; i < num_objects; ++i ) {
            for( int j = 0; j < num_objects; ++j ) {
                if( sd[i * num_objects + j] )
                    adj[i].emplace_back(j);
            }
        }

        // apply dfs starting from each vertex
        for( int r = 0; r < num_objects; ++r ) {
            std::vector<bool> visited(num_objects, false);
            std::vector<int> q = adj[r];
            while( !q.empty() ) {
                int i = q.back();
                q.pop_back();
                sd[r * num_objects + i] = true;
                if( !visited[i] ) {
                    visited[i] = true;
                    q.insert(q.end(), adj[i].begin(), adj[i].end());
                }
            }
        }
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            if( r == nullptr ) return nullptr; // CHECK: BLAI HACK
            assert((r != nullptr) && (r->size() == sample.num_states()));

            sample_denotation_t nr;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sr = (*r)[i];
                assert((sr != nullptr) && (sr->size() == state.num_objects() * state.num_objects()));

                state_denotation_t nsr(*sr);
                transitive_closure(state.num_objects(), nsr);
                nr.emplace_back(cache.find_or_insert_state_denotation(nsr));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nr, as_str()) : new sample_denotation_t(nr);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: PlusRole::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        // ATM let us call these Star(X) to get the same output than the Python module
        return std::string("Star(") + role_->as_str() + ")";
    }
};

class StarRole : public Role {
protected:
    const Role *role_;
    const PlusRole *plus_role_;

public:
    explicit StarRole(const Role *role)
            : Role(1 + role->complexity()),
              role_(role),
              plus_role_(new PlusRole(role)) {
    }
    ~StarRole() override {
        delete plus_role_;
    }
    [[nodiscard]] const Role* clone() const override {
        return new StarRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *pr = cache.find_sample_denotation(plus_role_->as_str());
            if( pr == nullptr ) return nullptr; // CHECK: BLAI HACK
            assert((pr != nullptr) && (pr->size() == sample.num_states()));

            sample_denotation_t nr;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sr = (*pr)[i];
                assert((sr != nullptr) && (sr->size() == state.num_objects() * state.num_objects()));

                state_denotation_t nsr(*sr);
                for( int j = 0; j < state.num_objects(); ++j ) {
                    int index = j * state.num_objects() + j;
                    nsr[index] = true;
                }
                nr.emplace_back(cache.find_or_insert_state_denotation(nsr));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nr, as_str()) : new sample_denotation_t(nr);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: StarRole::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] const Role* role() const { return role_; }


    [[nodiscard]] std::string as_str() const override {
        return std::string("Star(") + role_->as_str() + ")";
    }
};

class InverseRole : public Role {
protected:
    const Role *role_;

public:
    explicit InverseRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    ~InverseRole() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new InverseRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states()));

            sample_denotation_t nr;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sr = (*r)[i];
                assert((sr != nullptr) && (sr->size() == state.num_objects() * state.num_objects()));

                state_denotation_t nsr(state.num_objects() * state.num_objects(), false);
                for( int j = 0; j < state.num_objects(); ++j ) {
                    for( int k = 0; k < state.num_objects(); ++k ) {
                        int index = j * state.num_objects() + k;
                        if( (*sr)[index] ) {
                            int inv_index = k * state.num_objects() + j;
                            nsr[inv_index] = true;
                        }
                    }
                }
                nr.emplace_back(cache.find_or_insert_state_denotation(nsr));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nr, as_str()) : new sample_denotation_t(nr);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: InverseRole::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Inverse(") + role_->as_str() + ")";
    }
};

// RoleRestriction are only used for distance features
// and thus they are generated when generating such features
class RoleRestriction : public Role {
protected:
    const Role *role_;
    const Concept *restriction_;

public:
    RoleRestriction(const Role *role, const Concept *restriction)
            : Role(1 + role->complexity() + restriction->complexity()),
              role_(role),
              restriction_(restriction) {
    }
    ~RoleRestriction() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new RoleRestriction(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states()));
            const sample_denotation_t *d = cache.find_sample_denotation(restriction_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));

            sample_denotation_t nr;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *sr = (*r)[i];
                assert((sr != nullptr) && (sr->size() == state.num_objects() * state.num_objects()));
                const state_denotation_t *sd = (*d)[i];
                assert((sd != nullptr) && (sd->size() == state.num_objects()));

                state_denotation_t nsr(*sr);
                for( int j = 0; j < state.num_objects() * state.num_objects(); ++j ) {
                    if( nsr[j] ) {
                        //int src = j / state.num_objects();
                        int dst = j % state.num_objects();
                        nsr[j] = (*sd)[dst];
                    }
                }
                nr.emplace_back(cache.find_or_insert_state_denotation(nsr));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nr, as_str()) : new sample_denotation_t(nr);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: RoleRestriction::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Restrict(") + role_->as_str() + "," + restriction_->as_str() + ")";
    }
};

//! R - R' represents the set of pairs (a,b) such that R(a,b) and not R'(a,b)
//! Makes sense e.g. when R is a goal predicate
class RoleDifference : public Role {
protected:
    const Role *r1_;
    const Role *r2_;

public:
    RoleDifference(const Role *r1, const Role *r2)
            : Role(1 + r1->complexity() + r2->complexity()),
              r1_(r1),
              r2_(r2) {
    }
    ~RoleDifference() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new RoleDifference(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const override {
        const sample_denotation_t *cached = use_cache ? cache.find_sample_denotation(as_str()) : nullptr;
        if( cached == nullptr ) {
            const sample_denotation_t *r1_den = cache.find_sample_denotation(r1_->as_str());
            assert(r1_den && (r1_den->size() == sample.num_states()));
            const sample_denotation_t *r2_den = cache.find_sample_denotation(r2_->as_str());
            assert(r2_den && (r2_den->size() == sample.num_states()));

            sample_denotation_t nd;
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                unsigned m = state.num_objects();
                assert(state.id() == i);
                const state_denotation_t *sd1 = (*r1_den)[i];
                assert((sd1 != nullptr) && (sd1->size() == m*m));
                const state_denotation_t *sd2 = (*r2_den)[i];
                assert((sd2 != nullptr) && (sd2->size() == m*m));

                state_denotation_t nsd(m*m, false);

                for(int x = 0; x < m*m; ++x) {
                    if ((*sd1)[x] && !(*sd2)[x]) {
                        nsd[x] = true;
                    }
                }

                nd.emplace_back(cache.find_or_insert_state_denotation(nsd));
            }
            return use_cache ? cache.find_or_insert_sample_denotation(nd, as_str()) : new sample_denotation_t(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        throw std::runtime_error("Unexpected call: RoleDifference::denotation(Cache&, const Sample&, const State&)");
        return nullptr;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("RoleDifference(") + r1_->as_str() + "," + r2_->as_str() + ")";
    }
};

class Feature {
public:
    Feature() = default;
    virtual ~Feature() = default;
    virtual const Feature* clone() const = 0;

    [[nodiscard]] virtual int complexity() const = 0;
    [[nodiscard]] virtual int value(const Cache &cache, const Sample &sample, const State &state) const = 0;
    [[nodiscard]] virtual std::string as_str() const = 0;

    [[nodiscard]] std::string as_str_with_complexity() const {
        return std::to_string(complexity()) + "." + as_str();
    }

    friend std::ostream& operator<<(std::ostream &os, const Feature &f) {
        return os << f.as_str() << std::flush;
    }

    [[nodiscard]] virtual bool is_boolean() const = 0;
};

class NullaryAtomFeature : public Feature {
protected:
    const Predicate* predicate_;

public:
    explicit NullaryAtomFeature(const Predicate* predicate) : predicate_(predicate) { }
    ~NullaryAtomFeature() override = default;

    [[nodiscard]] const Feature* clone() const override {
        return new NullaryAtomFeature(*this);
    }

    [[nodiscard]] int complexity() const override { // Nullary atoms have complexity 0 by definition
        return 1;
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // Return true (1) iff some atom in the state has the coindiding predicate ID with the predicate of the
        // nullary atom, since there can only be one single atom with this predicate
        // TODO: This is unnecessarily expensive, and it is not cached anywhere
        for (const auto& atom_id:state.atoms()) {
            const Atom &atom = state.atom(atom_id);
            if( atom.is_instance(*predicate_)) {
                return 1;
            }
        }
        return 0;
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Atom[") + predicate_->name_ + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

class BooleanFeature : public Feature {
protected:
    const Concept *concept_;

public:
    explicit BooleanFeature(const Concept *concept) : Feature(), concept_(concept) { }
    ~BooleanFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new BooleanFeature(concept_);
    }

    [[nodiscard]] int complexity() const override {
        return concept_->complexity();
    }
    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // we look into cache for sample denotation using the concept name,
        // then index sample denotation with state id to find state denotation,
        // for finally computing cardinality (this assumes that state id is
        // index of state into sample.states())
        assert(sample.state(state.id()).id() == state.id());
        const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
        assert((d != nullptr) && (d->size() == sample.num_states()));
        const state_denotation_t *sd = (*d)[state.id()];
        assert((sd != nullptr) && (sd->size() == state.num_objects()));
        assert(sd->cardinality() < 2);
        return sd->cardinality();
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Bool[") + concept_->as_str() + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

class NumericalFeature : public Feature {
protected:
    const Concept *concept_;

public:
    explicit NumericalFeature(const Concept *concept) : Feature(), concept_(concept) { }
    ~NumericalFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new NumericalFeature(concept_);
    }

    [[nodiscard]] int complexity() const override {
        return concept_->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // we look into cache for sample denotation using the concept name,
        // then index sample denotation with state id to find state denotation,
        // for finally computing cardinality (this assumes that state id is
        // index of state into sample.states())
        assert(sample.state(state.id()).id() == state.id());
        const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
        assert((d != nullptr) && (d->size() == sample.num_states()));
        const state_denotation_t *sd = (*d)[state.id()];
        assert((sd != nullptr) && (sd->size() == state.num_objects()));
        return sd->cardinality();
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Num[") + concept_->as_str() + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return false; }
};

class DistanceFeature : public Feature {
protected:
    const Concept *start_;
    const Concept *end_;
    const Role *role_;

    mutable bool valid_cache_;
    mutable std::vector<int> cached_distances_;

public:
    DistanceFeature(const Concept *start, const Concept *end, const Role *role)
            : Feature(),
              start_(start),
              end_(end),
              role_(role),
              valid_cache_(false)
    {}
    ~DistanceFeature() override = default;
    const Feature* clone() const override {
        auto *f = new DistanceFeature(start_, end_, role_);
        f->valid_cache_ = valid_cache_;
        f->cached_distances_ = cached_distances_;
        return f;
    }

    int complexity() const override {
        return 1 + start_->complexity() + end_->complexity() + role_->complexity();
    }
    int value(const Cache &cache, const Sample &sample, const State &state) const override {
        assert(sample.state(state.id()).id() == state.id());
        const auto m = sample.num_states();
        if( !valid_cache_ ) {
            const sample_denotation_t *start_d = cache.find_sample_denotation(start_->as_str());
            assert((start_d != nullptr) && (start_d->size() == m));
            const sample_denotation_t *end_d = cache.find_sample_denotation(end_->as_str());
            assert((end_d != nullptr) && (end_d->size() == m));
            const sample_denotation_t *role_d = cache.find_sample_denotation(role_->as_str());
            assert((role_d != nullptr) && (role_d->size() == m));

            cached_distances_ = std::vector<int>(m, std::numeric_limits<int>::max());
            for( int i = 0; i < m; ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *start_sd = (*start_d)[i];
                assert((start_sd != nullptr) && (start_sd->size() == state.num_objects()));
                const state_denotation_t *end_sd = (*end_d)[i];
                assert((end_sd != nullptr) && (end_sd->size() == state.num_objects()));
                const state_denotation_t *role_sd = (*role_d)[i];
                assert((role_sd != nullptr) && (role_sd->size() == state.num_objects() * state.num_objects()));
                int distance = compute_distance(state.num_objects(), *start_sd, *end_sd, *role_sd);
                cached_distances_[i] = distance;
            }
            valid_cache_ = true;
        }
        return cached_distances_[state.id()];
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Dist[") + start_->as_str() + ";" + role_->as_str() + ";" + end_->as_str() + "]";
    }

    bool is_boolean() const override { return false; }
};


class ConditionalFeature : public Feature {
protected:
    const Feature* condition_;
    const Feature* body_;

public:
    ConditionalFeature(const Feature* condition, const Feature* body) :
            Feature(), condition_(condition), body_(body) { }

    ~ConditionalFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new ConditionalFeature(condition_, body_);
    }

    [[nodiscard]] int complexity() const override {
        return 1 + condition_->complexity() + body_->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        return value_from_components(cache, sample, state, condition_, body_);
    }

    static int value_from_components(const Cache &cache, const Sample &sample, const State &state,
                                     const Feature* condition, const Feature* body) {
        const auto condval = condition->value(cache, sample, state);
        return (condval) ? body->value(cache, sample, state) : std::numeric_limits<int>::max();
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("If{") + condition_->as_str() + "}{" + body_->as_str() + "}{Infty}";
    }

    [[nodiscard]] bool is_boolean() const override { return false; }
};


//! A feature with value f1 < f2
class DifferenceFeature : public Feature {
protected:
    const Feature* f1;
    const Feature* f2;

public:
    DifferenceFeature(const Feature* f1, const Feature* f2) :
            Feature(), f1(f1), f2(f2)
    {}

    ~DifferenceFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new DifferenceFeature(f1, f2);
    }

    [[nodiscard]] int complexity() const override {
        return 1 + f1->complexity() + f2->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        const auto f1val = f1->value(cache, sample, state);
        const auto f2val = f2->value(cache, sample, state);
        return f1val < f2val;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("LessThan{") + f1->as_str() + "}{" + f2->as_str() + "}";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

class Factory {
protected:
    const std::vector<std::string> nominals_;
    std::vector<const Role*> basis_roles_;
    std::vector<const Concept*> basis_concepts_;
    
    Options options;

    mutable std::vector<const Role*> roles_;

    // A layered set of concepts, concepts_[k] contains all concepts generated in the k-th application of
    // the concept grammar
    mutable std::vector<std::vector<const Concept*> > concepts_;
    std::vector<const Feature*> features_;

    //! Indices of features generated from goal concepts
    std::unordered_set<unsigned> goal_features_;

public:
    Factory(std::vector<std::string>  nominals, Options options) :
            nominals_(std::move(nominals)),
            options(std::move(options))
    {}
    virtual ~Factory() = default;

    void insert_basis(const Role *role) {
        basis_roles_.push_back(role);
    }
    void insert_basis(const Concept *concept) {
//        std::cout << "insert_basis: " << *concept << std::endl;
        basis_concepts_.push_back(concept);
    }

    void reset(bool remove) {
        // roles can be safely deleted since the
        // ones coming from the basis are clones
        while( !roles_.empty() ) {
            if( remove ) delete roles_.back();
            roles_.pop_back();
        }

        // concepts in all layers can be safely deleted since
        // the ones coming from the basis are clones
        while( !concepts_.empty() ) {
            if( remove ) {
                while( !concepts_.back().empty() ) {
                    delete concepts_.back().back();
                    concepts_.back().pop_back();
                }
            }
            concepts_.pop_back();
        }
        concepts_.pop_back();
    }

    static void insert_new_denotation_by_name(Cache &cache, const std::string &name, const sample_denotation_t *d) {
        cache.find_or_insert_sample_denotation_by_name(name, *d);
    }

    int insert_new_role(Cache &cache, const Role *role, const sample_denotation_t *d) const {
        roles_.push_back(role);
        cache.find_or_insert_sample_denotation(*d, role->as_str());
        return roles_.size();
    }

    //! Insert the given role as long as it is not redundant with some previous role and it is below the complexity
    //! bound. Return whether the role was effectively inserted.
    bool insert_role_if_possible(const Role &role, Cache &cache, const Sample &sample) const {
        auto sup = is_superfluous_or_exceeds_complexity_bound(role, cache, sample);
        if(!sup.first) {
//            std::cout << "INSERT: " + role.as_str() << std::endl;
            insert_new_role(cache, role.clone(), sup.second);
        } else {
            // We cannot just obviate this role, since some other role denotation may depend on this one's
//            std::cout << "PRUNE: " + role.as_str() << std::endl;
            insert_new_denotation_by_name(cache, role.as_str(), sup.second);
        }
        delete sup.second;
        return !sup.first;
    }

    static bool is_superfluous(Cache &cache, const sample_denotation_t *d) {
        return cache.find_sample_denotation(*d) != nullptr;
    }

    std::pair<bool, const sample_denotation_t*>
    is_superfluous_or_exceeds_complexity_bound(const Base &base, Cache &cache, const Sample &sample) const {
        if(base.complexity() > options.complexity_bound) {
//            std::cout << base.as_str() << " superfluous because complexity " << base.complexity() << ">" << options.complexity_bound << std::endl;
            return {true, nullptr};
        }

        const sample_denotation_t *d = base.denotation(cache, sample, false);
        const auto& denotation_idx = cache.cache1();
        auto it = denotation_idx.find(d);

        if (it != denotation_idx.end()) {
//            std::cout << base.as_str() << " (k=" << base.complexity() << ") superfluous because equivalent to " << it->second << std::endl;
            return {true, d};
        }
        return {false, d};
    }

    //!
    void attempt_concept_insertion(const Concept& concept, Cache& cache, const Sample& sample, int& pruning_count) const {
        if (concept.complexity() > options.complexity_bound) {
//            std::cout << concept.as_str() << " superfluous because complexity " << concept.complexity() << ">" << options.complexity_bound << std::endl;
            pruning_count++;
            return;
        }

        const sample_denotation_t *d = concept.denotation(cache, sample, false);
        const auto& index = cache.cache1();

        auto it = index.find(d);
        if (it != index.end()) {
            // Some other concept has the exact same denotation over the whole sample.
//            std::cout << concept.as_str() << " (k=" << concept.complexity() << ") superfluous because equivalent to " << it->second << std::endl;
            pruning_count++;

            // TODO Here we should be making sure that we're not pruning one simpler concept in favor of a
            //      more complex one. This is not trivial to do now due to the current architecture, because
            //      if we want to replace an already-existing concept with a new one, we should guarantee that there
            //      are no third concepts that have been generated in the meantime and make use of the one we want
            //      to eliminate

        } else {
            assert(!concepts_.empty());
            concepts_.back().push_back(concept.clone());
            cache.find_or_insert_sample_denotation(*d, concept.as_str());
        }

        delete d;
    }

    bool check_timeout(const std::clock_t& start_time) const {
        if (options.timeout <= 0) return false;
        double elapsed = (std::clock() - start_time) / (double) CLOCKS_PER_SEC;
        if (elapsed > options.timeout) {
            std::cout << "\tTimeout of " << options.timeout << " sec. reached while generating concepts" << std::endl;
            return true;
        }
        return false;
    }

    // apply one iteration of the concept generation grammar
    // new concepts are left on a new last layer in concepts_
    // new concepts are non-redundant if sample != nullptr
    int advance_step(Cache &cache, const Sample &sample, const std::clock_t& start_time) const {
        int num_pruned_concepts = 0;
        if( concepts_.empty() ) { // On the first iteration, we simply process the basis concepts and return
            concepts_.emplace_back();
            for (const auto* concept : basis_concepts_) {
                attempt_concept_insertion(*concept, cache, sample, num_pruned_concepts);
            }
            return num_pruned_concepts;
        }


        // Classify concepts and roles by their complexity - we will use this to minimize complexity checks later
        std::vector<std::vector<const Role*>> roles_by_complexity(options.complexity_bound+1);
        std::vector<std::vector<const Concept*>> concepts_in_last_layer_by_complexity(options.complexity_bound+1);
        std::vector<std::vector<const Concept*>> concepts_in_previous_layers_by_complexity(options.complexity_bound+1);

        for (const auto* r:roles_) roles_by_complexity.at(r->complexity()).push_back(r);  // TODO Do this just once
        for (const auto* c:concepts_.back()) concepts_in_last_layer_by_complexity.at(c->complexity()).push_back(c);

        for( unsigned layer = 0; layer < concepts_.size()-1; ++layer ) {
            for (const auto* c:concepts_[layer]) {
                concepts_in_previous_layers_by_complexity.at(c->complexity()).push_back(c);
            }
        }

        // create new concept layer
        bool is_first_non_basis_iteration = (concepts_.size() == 1);
        concepts_.emplace_back();

        if (is_first_non_basis_iteration) {
            // Let's create equal concepts for those pairs of roles (whose number is already fixed at this point)
            // such that both arise from same predicate and one of them is the "goal version" of the other, e.g.
            // on and on_g, in blocksworld.
            for (const auto* r1:roles_) {
                for (const auto* r2:roles_) {
                    if (dynamic_cast<const RoleDifference*>(r1) || dynamic_cast<const RoleDifference*>(r2)) {
                        // R=R' Makes no sense when R or R' are already a role difference
                        continue;
                    }
                    auto name1 = get_role_predicate(r1)->name();
                    auto name2 = get_role_predicate(r2)->name();

                    // A hacky way to allow only equal concepts R=R' where R and R' come from the same predicate
                    // and the first one is a goal role
                    if (name1.substr(name1.size()-2) != "_g") continue;
                    if (name1.substr(0, name1.size()-2) != name2) continue;

                    EqualConcept eq_concept(r1, r2);
                    // Force the complexity of R=R' to be 1, if both roles are of the same type
                    // This is currently disabled, as the concept-pruning strategy doesn't work as expected when
                    // this is enabled - sometimes concepts with lower complexity get pruned in favor of others
                    // with higher complexity
                    // if (typeid(*r1) == typeid(*r2)) eq_concept.force_complexity(1);

                    attempt_concept_insertion(eq_concept, cache, sample, num_pruned_concepts);

                    if (check_timeout(start_time)) return -1;
                }
            }
        }

        for (int k = 0; k <= (options.complexity_bound-1); ++k) {
            for (const auto* concept:concepts_in_last_layer_by_complexity[k]) {

                // Negate concepts in the last layer, only if they are not already negations
                if (!dynamic_cast<const NotConcept*>(concept)) {
                    attempt_concept_insertion(NotConcept(concept), cache, sample, num_pruned_concepts);
                }


                // generate exist and forall combining a role with a concept in the last layer
                for (int k2 = 0; k2 <= (options.complexity_bound-k-1); ++k2) {
                    for (const auto *role:roles_by_complexity[k2]) {
                        attempt_concept_insertion(ExistsConcept(concept, role), cache, sample, num_pruned_concepts);
                        attempt_concept_insertion(ForallConcept(concept, role), cache, sample, num_pruned_concepts);

                        if (check_timeout(start_time)) return -1;
                    }
                }
                if (check_timeout(start_time)) return -1;
            }
        }


        // generate conjunctions of (a) a concept of the last layer with a concept in some previous layer, and
        // (b) two concepts of the last layer, avoiding symmetries
        for (int k = 0; k <= (options.complexity_bound-1); ++k) {
            for (unsigned i_k = 0; i_k < concepts_in_last_layer_by_complexity[k].size(); ++i_k) {
                const auto& concept1 = concepts_in_last_layer_by_complexity[k][i_k];

                // (a) a concept of the last layer with a concept in some previous layer
                for (int k2 = 0; k2 <= (options.complexity_bound-k); ++k2) {
                    for (const auto *concept2:concepts_in_previous_layers_by_complexity[k2]) {
                        attempt_concept_insertion(AndConcept(concept1, concept2), cache, sample, num_pruned_concepts);

                        if (check_timeout(start_time)) return -1;
                    }
                }


                // (b) two concepts of the last layer, avoiding symmetries
                for (int k2 = k; k2 <= (options.complexity_bound-k); ++k2) {
                    unsigned start_at = (k == k2) ? i_k+1: 0; // Break symmetries within same complexity bucket
                    for (unsigned i_k2 = start_at; i_k2 < concepts_in_last_layer_by_complexity[k2].size(); ++i_k2) {
                        const auto& concept2 = concepts_in_last_layer_by_complexity[k2][i_k2];
                        attempt_concept_insertion(AndConcept(concept1, concept2), cache, sample, num_pruned_concepts);

                        if (check_timeout(start_time)) return -1;
                    }
                }
            }
        }
        return num_pruned_concepts;
    }

    //! Retrieve the predicate at the basis of a given role (given the current grammar restrictions, there will be
    //! exactly one such predicate
    static const Predicate* get_role_predicate(const Role* r) {

        if (const auto *c = dynamic_cast<const PrimitiveRole*>(r)) {
            return c->predicate();

        } else if (const auto* c = dynamic_cast<const StarRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const PlusRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const RoleRestriction*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const InverseRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const RoleDifference*>(r)) {
            throw std::runtime_error("Unimplemented");
        }
        throw std::runtime_error("Unknown role type");
    }

    void generate_basis(const Sample &sample) {
        insert_basis(new UniversalConcept);
        insert_basis(new EmptyConcept);
        for( int i = 0; i < int(sample.num_predicates()); ++i ) {
            const Predicate *predicate = &sample.predicates()[i];
            if( predicate->arity() == 1 ) {
                insert_basis(new PrimitiveConcept(predicate));
            } else if( predicate->arity() == 2 ) {
                insert_basis(new PrimitiveRole(predicate));
            }
        }

        for (const auto& nominal:nominals_) {
            insert_basis(new NominalConcept(nominal));
        }
        std::cout << "BASIS: #concepts=" << basis_concepts_.size() << ", #roles=" << basis_roles_.size() << std::endl;
//        report_dl_data(std::cout);
    }

    std::vector<const Concept*> generate_goal_concepts_and_roles(Cache &cache, const Sample &sample) {
        std::vector<const Concept*> goal_concepts;

        const Sample::PredicateIndex& predicate_idx = sample.predicate_index();
        for (const predicate_id_t& pid:sample.goal_predicates()) {
            const Predicate& pred = sample.predicate(pid);
            unsigned arity = pred.arity();
            const auto& name = pred.name();
            auto name_n = name.size();
            assert(name[name_n-2] == '_' && name[name_n-1] == 'g');  // Could try something more robust, but ATM this is OK
            auto nongoal_name = name.substr(0, name_n-2);
            const Predicate& nongoal_pred = sample.predicate(predicate_idx.at(nongoal_name));
            if (arity == 0) {
                // No need to do anything here, as the corresponding NullaryAtomFeature will always be generated
                // with complexity 1
                assert (0); // This will need some refactoring
//                goal_features.push_back(new NullaryAtomFeature(&pred));

            } else if (arity == 1) {
                // Force into the basis a concept "p_g AND Not p"
                auto *c = new PrimitiveConcept(&nongoal_pred);
                auto *not_c = new NotConcept(c);
                auto *c_g = new PrimitiveConcept(&pred);
                auto *and_c = new AndConcept(c_g, not_c);

                // Insert the denotations into the cache
//                for (const auto elem:std::vector<const Concept*>({and_c})) {
                for (const auto* elem:std::vector<const Concept*>({c, not_c, c_g, and_c})) {
                    const sample_denotation_t *d = elem->denotation(cache, sample, false);
                    cache.find_or_insert_sample_denotation(*d, elem->as_str());
                }
//                for (const auto elem:std::vector<const Concept*>({c, not_c, c_g})) {
//                    cache.remove_sample_denotation(elem->as_str());
//                }


                and_c->force_complexity(1);
                for (const auto* elem:std::vector<const Concept*>({c, not_c, c_g, and_c})) {
                    insert_basis(elem);
                }

//                goal_features.push_back(new NumericalFeature(and_c));
                goal_concepts.push_back(and_c);

            } else if (arity == 2) {
                // Force into the basis a new RoleDifference(r_g, r)
                auto *r = new PrimitiveRole(&nongoal_pred);
                auto *r_g = new PrimitiveRole(&pred);
                auto *r_diff = new RoleDifference(r_g, r);

                // Insert the denotations into the cache
//                for (const auto elem:std::vector<const Role*>({r_diff})) {
//                for (const auto elem:std::vector<const Role*>({r, r_g, r_diff})) {
//                    const sample_denotation_t *d = elem->denotation(cache, sample, false);
//                    cache.find_or_insert_sample_denotation(*d, elem->as_str());
//                }
//                for (const auto elem:std::vector<const Role*>({r, r_g})) {
//                    cache.remove_sample_denotation(elem->as_str());
//                }

                r_diff->force_complexity(1);
                for (const auto* elem:std::vector<const Role*>({r, r_g, r_diff})) {
                    insert_basis(elem);
                }

                auto *ex_c = new ExistsConcept(new UniversalConcept, r_diff); // memleak
//                goal_features.push_back(new NumericalFeature(ex_c));
                goal_concepts.push_back(ex_c);
            }
        }
        return goal_concepts;
    }

    int generate_roles(Cache &cache, const Sample &sample) const {
        assert(roles_.empty());

        // Insert the basis (i.e. primitive) roles as long as they are not redundant
        for (const auto *role : basis_roles_) {
            auto p = is_superfluous_or_exceeds_complexity_bound(*role, cache, sample);
            if( !p.first ) {
                insert_new_role(cache, role->clone(), p.second);
            } else {
                // we cannot just obviate this role since another
                // role denotation may depend on this denotation
                //std::cout << "PRUNE: " + role.as_str() << std::endl;
                insert_new_denotation_by_name(cache, role->as_str(), p.second);
            }
            delete p.second;
        }

        for (const auto *role : basis_roles_) {
            // Create Plus(R) roles from the primitive roles
            PlusRole p_role(role);
            insert_role_if_possible(p_role, cache, sample);

            // Create Star(R) roles from the primitive roles !!! NOTE ATM we deactivate Star roles
            // insert_role_if_possible(StarRole(role), cache, sample);

            // Create Inverse(R) roles from the primitive roles
            InverseRole inv_r(role);
            insert_role_if_possible(InverseRole(role), cache, sample);

            // Create Inverse(Plus(R)) and Plus(Inverse(R)) roles from the primitive roles
            insert_role_if_possible(InverseRole(p_role.clone()), cache, sample);
            insert_role_if_possible(PlusRole(inv_r.clone()), cache, sample);
        }
        std::cout << "ROLES: #roles=" << roles_.size() << std::endl;
        return roles_.size();
    }

    std::vector<const Concept*> generate_concepts(Cache &cache, const Sample &sample, const std::clock_t& start_time) const {
        std::size_t num_concepts = 0;
        bool some_new_concepts = true;
        bool timeout_reached = false;
        for( int iteration = 0; some_new_concepts && !timeout_reached; ++iteration ) {
            std::cout << "DL::concept-generation: iteration=" << iteration
                      << ", #concepts=" << num_concepts
                      << ", #concepts-in-last-layer=" << (concepts_.empty() ? 0 : concepts_.back().size())
                      << std::endl;

            std::size_t num_concepts_before_step = num_concepts;

            // Let the fun begin:
            auto num_pruned_concepts = advance_step(cache, sample, start_time);
            timeout_reached = num_pruned_concepts < 0;

            num_concepts += concepts_.empty() ? 0 : concepts_.back().size();
            some_new_concepts = num_concepts > num_concepts_before_step;

            std::cout << "\tResult of iteration:"
                      << " #concepts-in-layer=" << concepts_.back().size()
                      << ", #pruned-concepts=" << num_pruned_concepts
                      << std::endl;
//            report_dl_data(std::cout);
        }

        if (concepts_.back().empty()) {  // i.e. we stopped because last iteration was already a fixpoint
            concepts_.pop_back();
        }

        // Important to use this vector from here on, as it is sorted by complexity
        auto all_concepts = consolidate_concepts();
        assert(all_concepts.size() == num_concepts);

        std::cout << "DL::Factory: #concepts-final=" << num_concepts << std::endl;
        return all_concepts;
    }

    void generate_comparison_features(
            const std::vector<const Feature*>& base_features,
            Cache& cache,
            const Sample& sample,
            feature_cache_t& seen_denotations)
    {
        if (!options.comparison_features) return;

        // get the max index here, as we'll keep adding elements to the same `features_` vector:
        auto n = features_.size();

        for (std::size_t i = 0; i < n; ++i) {
            const auto* f_i = features_[i];
            if (f_i->is_boolean() || f_i->complexity() + 1 + 1 > options.complexity_bound) continue;

            for (std::size_t j = 0; j < n; ++j) {
                const auto* f_j = features_[j];
                if (i == j || f_j->is_boolean() || f_i->complexity() + f_j->complexity() + 1 > options.complexity_bound)
                    continue;

                const auto *feature = new DifferenceFeature(f_i, f_j);
                if (!insert_feature_if_necessary(
                        feature, options.complexity_bound, cache, sample, seen_denotations)) {
                    delete feature;
                }
            }
        }
    }

    void generate_conditional_features(
            const std::vector<const Feature*>& base_features,
            Cache& cache,
            const Sample& sample,
            feature_cache_t& seen_denotations)
    {
        if (options.cond_complexity_bound <= 0) return;

        for (std::size_t i = 0, n = features_.size(); i < n; ++i) {
            const auto* cond = features_[i];
            if (!cond->is_boolean() || cond->complexity() + 1 + 1 > options.cond_complexity_bound) continue;

            for (std::size_t j = i + 1; j < n; ++j) {
                const auto* body = features_[j];
                if (body->is_boolean() ||
                    cond->complexity() + body->complexity() + 1 > options.cond_complexity_bound)
                    continue;

                const auto *feature = new ConditionalFeature(cond, body);
                if (!insert_feature_if_necessary(
                        feature, options.cond_complexity_bound, cache, sample, seen_denotations)) {
                    delete feature;
                }
            }
        }
    }

    void generate_features(
            const std::vector<const Concept*>& concepts,
            Cache &cache, const Sample &sample,
            const ::Sample::TransitionSample& transitions,
            const std::vector<const Concept*>& forced_goal_features);

    void print_feature_count() const {
        unsigned num_nullary_features = 0, num_boolean_features = 0, num_numeric_features = 0,
                num_distance_features = 0, num_conditional_features = 0, num_comparison_features = 0;
        auto nf = features_.size();
        for (const auto *f:features_) {
            if (dynamic_cast<const NullaryAtomFeature*>(f)) num_nullary_features++;
            else if (dynamic_cast<const BooleanFeature*>(f)) num_boolean_features++;
            else if (dynamic_cast<const NumericalFeature*>(f)) num_numeric_features++;
            else if (dynamic_cast<const DistanceFeature*>(f)) num_distance_features++;
            else if (dynamic_cast<const ConditionalFeature*>(f)) num_conditional_features++;
            else if (dynamic_cast<const DifferenceFeature*>(f)) num_comparison_features++;
            else throw std::runtime_error("Unknown feature type");
        }

        assert(nf == num_nullary_features+num_boolean_features+num_numeric_features
                     +num_distance_features+num_conditional_features+num_comparison_features);
        std::cout << "FEATURES: #features=" << nf
                  << ", #nullary="   << num_nullary_features
                  << ", #boolean="   << num_boolean_features
                  << ", #numerical=" << num_numeric_features
                  << ", #distance="  << num_distance_features
                  << ", #conditional="  << num_conditional_features
                  << ", #comparison="  << num_comparison_features
                  << std::endl;
    }

    bool generate_cardinality_feature_if_not_redundant(
            const Concept* c,
            Cache &cache,
            const Sample &sample,
            const ::Sample::TransitionSample& transitions,
            feature_cache_t& seen_denotations,
            bool can_be_pruned);

    //! Insert the given feature if its complexity is below the given bound, its denotation is not constant,
    //! and its denotation trail is not redundant with that of some previously-generated feature
    bool insert_feature_if_necessary(
            const Feature* feature, unsigned bound,
            Cache &cache, const Sample &sample, feature_cache_t& seen_denotations);

    void generate_distance_features(const std::vector<const Concept*>& concepts, Cache &cache, const Sample &sample, feature_cache_t& seen_denotations) {
        if (options.dist_complexity_bound<=0) return;

        auto m = sample.num_states();

        // identify concepts with singleton denotation across all states
        // as these are the candidates for start concepts
        std::vector<const Concept*> start_concepts;
        for (const Concept* c:concepts) {
            const sample_denotation_t *d = cache.find_sample_denotation(c->as_str());
            assert((d != nullptr) && (d->size() == m));
            bool singleton_denotations = true;
            for( int j = 0; singleton_denotations && (j < m); ++j )
                singleton_denotations = (*d)[j]->cardinality() == 1;
            if( singleton_denotations ) {
                start_concepts.push_back(c);
            }
        }

        // create role restrictions to be used in distance features
        std::vector<const Role*> dist_roles(roles_);  // Start with all base roles
        for (const Role* r:roles_) {
            for (const Concept* c:concepts) {
                RoleRestriction role_restriction(r, c);
                if( role_restriction.complexity()+3 > options.dist_complexity_bound ) continue;
                const sample_denotation_t *d = role_restriction.denotation(cache, sample, false);
                if( !is_superfluous(cache, d) ) {
                    assert(cache.cache1().find(d) == cache.cache1().end());
                    dist_roles.push_back(role_restriction.clone());
                    cache.find_or_insert_sample_denotation(*d, role_restriction.as_str());
                    assert(cache.cache1().find(d) != cache.cache1().end());
                    //std::cout << "ACCEPT RR(sz=" << cache_for_role_restrictions.cache1().size() << "): " + role_restriction.as_str_with_complexity() << std::endl;
                } else {
                    //std::cout << "PRUNE RR: " + role_restriction.as_str() << std::endl;
                }
                delete d;
            }
        }

        // create distance features
        int num_distance_features = 0;
        for(const Concept* start:start_concepts) {
            for (const Concept* end:concepts) {
                if (start == end) continue;
                for (const Role* role:dist_roles) {
                    DistanceFeature df(start, end, role);
                    if( df.complexity() > options.dist_complexity_bound ) continue;
                    //std::cout << "TESTING: " << df.as_str_with_complexity() << std::endl;

                    // fill cache with denotations for start and end concepts
                    const sample_denotation_t *ds = cache.find_sample_denotation(start->as_str());
                    assert((ds != nullptr) && (ds->size() == m));
                    cache.find_or_insert_sample_denotation(*ds, start->as_str());
                    const sample_denotation_t *de = cache.find_sample_denotation(end->as_str());
                    assert((de != nullptr) && (de->size() == m));
                    cache.find_or_insert_sample_denotation(*de, end->as_str());

                    // check whether df is novel
                    feature_sample_denotation_t fd(m);
                    for( int index = 0; index < m; ++index ) {
                        assert(sample.state(index).id() == index);
                        const State &state = sample.state(index);
                        int value = df.value(cache, sample, state);
                        fd[index] = value;
                    }

                    auto it = seen_denotations.find(fd);
                    if( it == seen_denotations.end() ) {
                        ++num_distance_features;
                        features_.emplace_back(df.clone());
                        seen_denotations.emplace(fd, features_.back());
                        //std::cout << "ACCEPT: " + df.as_str_with_complexity() << std::endl;
                    } else {
                        //std::cout << "REJECT: " + df.as_str() << std::endl;
                        //std::cout << "PRUNED-BY: " << it->second->as_str() << std::endl;
                    }
                }
            }
        }
    }

    std::ostream& report_dl_data(std::ostream &os) const {

        unsigned nconcepts = 0;
        for (const auto& layer:concepts_) nconcepts += layer.size();
        auto nroles = roles_.size();
        os << "Total number of concepts / roles: " << nconcepts << "/" << nroles << std::endl;

        os << "Base concepts (sz=" << basis_concepts_.size() << "): ";
        for( int i = 0; i < int(basis_concepts_.size()); ++i ) {
            os << basis_concepts_[i]->as_str();
            if( 1 + i < int(basis_concepts_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "Base roles (sz=" << basis_roles_.size() << "): ";
        for( int i = 0; i < int(basis_roles_.size()); ++i ) {
            os << basis_roles_[i]->as_str();
            if( 1 + i < int(basis_roles_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "All Roles under complexity " << options.complexity_bound << " (sz=" << roles_.size() << "): ";
        for( int i = 0; i < int(roles_.size()); ++i ) {
            os << roles_[i]->as_str_with_complexity();
            if( 1 + i < int(roles_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "All concepts (by layer) under complexity " << options.complexity_bound << ": " << std::endl;
        for( int layer = 0; layer < concepts_.size(); ++layer ) {
            os << "    Layer " << layer << " (sz=" << concepts_[layer].size() << "): ";
            for( int i = 0; i < int(concepts_[layer].size()); ++i ) {
                os << concepts_[layer][i]->as_str_with_complexity();
                if( 1 + i < int(concepts_[layer].size()) ) os << ", ";
            }
            os << std::endl;
        }

        os << "Nullary-atom features: ";
        bool need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const NullaryAtomFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Boolean features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const BooleanFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Numerical features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const NumericalFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Distance features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const DistanceFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        return os;
    }

    void output_feature_matrix(std::ostream &os, const Cache &cache, const Sample &sample) const {
        auto num_features = features_.size();
        for( int i = 0; i < sample.num_states(); ++i ) {
            const State &state = sample.state(i);

            // one line per state with the numeric denotation of all features
            for( int j = 0; j < num_features; ++j ) {
                const Feature &feature = *features_[j];
                os << feature.value(cache, sample, state);
                if( 1 + j < num_features ) os << " ";
            }
            os << std::endl;
        }
    }

    void output_feature_info(std::ostream &os, const Cache &cache, const Sample &sample) const {
        auto num_features = features_.size();

        // Line #1: feature names
        for( int i = 0; i < num_features; ++i ) {
            const Feature &feature = *features_[i];
            os << feature.as_str();
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #2: feature complexities
        for( int i = 0; i < num_features; ++i ) {
            const Feature &feature = *features_[i];
            os << feature.complexity();
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #3: feature types (0: boolean; 1: numeric)
        for( int i = 0; i < num_features; ++i ) {
            const Feature* feature = features_[i];
            os << (feature->is_boolean() ? 0 : 1);
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #4: whether feature is goal feature (0: No; 1: yes)
        for( int i = 0; i < num_features; ++i ) {
            auto it = goal_features_.find(i);
            os << (it == goal_features_.end() ? 0 : 1);
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;
    }

    void log_all_concepts_and_features(const std::vector<const Concept*>& concepts,
                                       const SLTP::DL::Cache &cache, const SLTP::DL::Sample &sample,
                                       const std::string& workspace, bool print_denotations);

    //! Return all generated concepts in a single, *unlayered* vector, and sorted by increasing complexity
    std::vector<const Concept*> consolidate_concepts() const {
        std::vector<const Concept*> all;
        for (const auto& layer:concepts_) all.insert(all.end(), layer.begin(), layer.end());

        std::sort(std::begin(all), std::end(all), [](const Concept* c1, const Concept* c2) {
            return c1->complexity() < c2->complexity();
        });

        return all;
    }
};

} // namespaces
