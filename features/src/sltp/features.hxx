
#ifndef FEATURES_HXX
#define FEATURES_HXX

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


#include <sltp/utils.hxx>

namespace SLTP {

namespace DL {

using object_id_t = unsigned;
using predicate_id_t = unsigned;
using atom_id_t = unsigned;
using state_id_t = unsigned;

class Atom;
class Sample;
class Concept;
class Role;
class State;

// A state denotation for concept C is a subset of objects
// implemented as a bitmap (ie. vector<bool>). These are
// cached so that at most one copy of a bitmap exists.
struct state_denotation_t : public std::vector<bool> {
    state_denotation_t(size_t n, bool value) : std::vector<bool>(n, value) { }
    size_t cardinality() const {
        size_t n = 0;
        for( size_t i = 0; i < size(); ++i )
            n += (*this)[i];
        return n;
    }
};

// A (full) sample denotation for concept C is a vector of
// state denotations, one per each state in the sample.
// Since we cache state denotations, a sample denotation
// is implemented as a vector of pointers to bitmaps.
using sample_denotation_t = std::vector<const state_denotation_t*>;

using feature_state_denotation_t = int;
using feature_sample_denotation_t = std::vector<feature_state_denotation_t>;

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
            for( int i = 0; i < int(obj->size()); ++i )
                hash = hash ^ (*this)((*obj)[i]);
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
    using cache4_t = std::unordered_map<std::string, const Atom*>;

  protected:
    cache1_t cache1_;
    cache2_t cache2_;
    cache3_t cache3_;
    cache4_t cache4_;

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
        cache1_t::const_iterator it = cache1_.find(&d);
        return it == cache1_.end() ? nullptr : it->first;
    }
    const sample_denotation_t* find_or_insert_sample_denotation(const sample_denotation_t &d, const std::string &name) {
        cache1_t::const_iterator it = cache1_.find(&d);
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
        cache2_t::const_iterator it = cache2_.find(name);
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
        cache2_t::const_iterator it = cache2_.find(name);
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
    const state_denotation_t* find_state_denotation(const state_denotation_t &sd) const {
        cache3_t::const_iterator it = cache3_.find(&sd);
        return it == cache3_.end() ? nullptr : *it;
    }

    const state_denotation_t* find_or_insert_state_denotation(const state_denotation_t &sd) {
        cache3_t::const_iterator it = cache3_.find(&sd);
        if( it == cache3_.end() ) {
            const state_denotation_t *nsd = new state_denotation_t(sd);
            cache3_.insert(nsd);
            return nsd;
        } else {
            return *it;
        }
    }

    // cache4: atoms
    //
    // name -> atom
    const Atom* find_atom(const std::string &name) const {
        cache4_t::const_iterator it = cache4_.find(name);
        return it == cache4_.end() ? nullptr : it->second;
    }
    const Atom* find_or_insert_atom(const Atom &atom);

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
    Object(unsigned id, const std::string &name) : id_(id), name_(name) { }
    int id() const {
        return id_;
    }
    const std::string& as_str() const {
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
    Predicate(unsigned id, const std::string &name, int arity)
      : id_(id), name_(name), arity_(arity) {
    }
    predicate_id_t id() const {
        return id_;
    }
    int arity() const {
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

    predicate_id_t pred_id() const {
        return predicate_;
    }
    const std::vector<object_id_t>& objects() const {
        return objects_;
    }

    // Return the i-th object of the current atom
    object_id_t object(int i) const {
        return objects_.at(i);
    }

    bool is_instance(const Predicate &predicate) const {
        return predicate_ == predicate.id();
    }

    std::vector<unsigned> data() const {
        std::vector<unsigned> res(1, predicate_);
        res.insert(res.end(), objects_.begin(), objects_.end());
        return res;
    }

    std::string as_str(const Sample &sample) const;
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

  protected:
    const std::string name_;
    const std::vector<Object> objects_;
    const std::vector<Atom> atoms_;

    // mapping from object names to their ID in the sample
    ObjectIndex object_index_;

    // mapping from <predicate name, obj_name, ..., obj_name> to the ID of the corresponding GroundPredicate
    AtomIndex atom_index_;

  public:
    Instance(const std::string &name,
             std::vector<Object> &&objects,
             std::vector<Atom> &&atoms,
             ObjectIndex &&object_index,
             AtomIndex &&atom_index)
      : name_(name),
        objects_(std::move(objects)),
        atoms_(std::move(atoms)),
        object_index_(std::move(object_index)),
        atom_index_(std::move(atom_index)) {
    }
    Instance(const Instance& ins) = default;
    Instance(Instance &&ins) = default;
    ~Instance() = default;

    const std::string& name() const {
        return name_;
    }
    int num_objects() const {
        return objects_.size();
    }
    int num_atoms() const {
        return atoms_.size();
    }
    const Atom& atom(unsigned id) const {
        return atoms_.at(id);
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

    unsigned id() const {
        return id_;
    }
    const std::vector<atom_id_t>& atoms() const {
        return atoms_;
    }
    int num_objects() const {
        return instance_.num_objects();
    }

    const Instance& instance() const {
        return instance_;
    }
    const Atom& atom(atom_id_t id) const {
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
    const std::string name_;
    const std::vector<Predicate> predicates_;
    const std::vector<Instance> instances_;
    const std::vector<State> states_;

    // mapping from predicate names to their ID in the sample
    PredicateIndex predicate_index_;

    Sample(std::string &name,
           std::vector<Predicate> &&predicates,
           std::vector<Instance> &&instances,
           std::vector<State> &&states,
           PredicateIndex &&predicate_index)
      : name_(name),
        predicates_(std::move(predicates)),
        instances_(std::move(instances)),
        states_(std::move(states)),
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

    const std::string& name() const {
        return name_;
    }
    int num_predicates() const {
        return predicates_.size();
    }
    int num_states() const {
        return states_.size();
    }

    const std::vector<Predicate>& predicates() const {
        return predicates_;
    }
    const std::vector<Instance>& instances() const {
        return instances_;
    }
    const std::vector<State>& states() const {
        return states_;
    }

    const Instance& instance(int i) const {
        return instances_.at(i);
    }
    const Predicate& predicate(unsigned id) const {
        return predicates_.at(id);
    }
    const State& state(unsigned id) const {
        return states_.at(id);
    }

    // factory method - reads sample from serialized data
    static const Sample read(std::istream &is);
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
    // unsigned id_;
    int complexity_;

  public:
    explicit Base(int complexity) : complexity_(complexity) { }
    // const unsigned id() const { return id_; }
    const int complexity() const { return complexity_; }

    //virtual void denotation(Denotation &d) const = 0;
    virtual const sample_denotation_t* denotation(Cache &cache, const Sample &sample, bool use_cache) const = 0;
    virtual const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const = 0;
    virtual std::string as_str() const = 0;

    std::string as_str_with_complexity() const {
        return std::to_string(complexity_) + "." + as_str();
    }

    friend std::ostream& operator<<(std::ostream &os, const Base &base) {
        return os << base.as_str_with_complexity() << std::flush;
    }
};

class Concept : public Base {
  public:
    explicit Concept(int complexity) : Base(complexity) { }
    virtual ~Concept() = default;
    virtual const Concept* clone() const = 0;
};

class Role : public Base {
  public:
    explicit Role(int complexity) : Base(complexity) { }
    virtual ~Role() = default;
    virtual const Role* clone() const = 0;
};

class PrimitiveConcept : public Concept {
  protected:
    const Predicate *predicate_;

  public:
    explicit PrimitiveConcept(const Predicate *predicate) : Concept(1), predicate_(predicate) { }
    ~PrimitiveConcept() override = default;
    const Concept* clone() const override {
        return new PrimitiveConcept(predicate_);
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
    std::string as_str() const override {
        return predicate_->name_;
    }
};


class NominalConcept : public Concept {
protected:
    const std::string& name_;

public:
    explicit NominalConcept(const std::string &name) : Concept(1), name_(name) {}
    ~NominalConcept() override = default;

    const Concept* clone() const override {
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
    std::string as_str() const override {
        return std::string("Nominal(") + name_ + ")";
    }
};

class UniversalConcept : public Concept {
public:
    UniversalConcept() : Concept(0) { }
    ~UniversalConcept() override = default;
    const Concept* clone() const override {
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
    std::string as_str() const override {
        return "<universe>";
    }
};


class EmptyConcept : public Concept {
  public:
    EmptyConcept() : Concept(0) { }
    ~EmptyConcept() override = default;
    const Concept* clone() const override {
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
    std::string as_str() const override {
        return "<empty>";
    }
};

class AndConcept : public Concept {
  protected:
    const Concept *concept1_;
    const Concept *concept2_;

  public:
    AndConcept(const Concept *concept1, const Concept *concept2)
      : Concept(1 + concept1->complexity() + concept2->complexity()),
        concept1_(concept1),
        concept2_(concept2) {
    }
    ~AndConcept() override = default;
    const Concept* clone() const override {
        return new AndConcept(concept1_, concept2_);
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

    std::string as_str() const override {
        return std::string("And(") +  concept1_->as_str() + "," + concept2_->as_str() + ")";
    }
};

class NotConcept : public Concept {
  protected:
    const Concept *concept_;

  public:
    explicit NotConcept(const Concept *concept)
      : Concept(1 + concept->complexity()),
        concept_(concept) {
    }
    ~NotConcept() override = default;
    const Concept* clone() const override {
        return new NotConcept(concept_);
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

    std::string as_str() const override {
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
    const Concept* clone() const override {
        return new ExistsConcept(concept_, role_);
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

    std::string as_str() const override {
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
    const Concept* clone() const override {
        return new ForallConcept(concept_, role_);
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

    std::string as_str() const override {
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
    const Concept* clone() const override { return new EqualConcept(*this); }

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
    std::string as_str() const override {
        return std::string("Equal(") + r1_->as_str() + "," + r2_->as_str() + ")";
    }
};


class PrimitiveRole : public Role {
  protected:
    const Predicate *predicate_;

  public:
    explicit PrimitiveRole(const Predicate *predicate) : Role(1), predicate_(predicate) { }
    ~PrimitiveRole() override { }
    const Role* clone() const override {
        return new PrimitiveRole(predicate_);
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

    std::string as_str() const override {
        return predicate_->name_;
    }
};

class PlusRole : public Role {
  protected:
    const Role *role_;

  public:
    explicit PlusRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    ~PlusRole() override { }
    const Role* clone() const override {
        return new PlusRole(role_);
    }

    // apply Johnson's algorithm for transitive closure
    void transitive_closure(int num_objects, state_denotation_t &sd) const {
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

    std::string as_str() const override {
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
    const Role* clone() const override {
        return new StarRole(role_);
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

    std::string as_str() const override {
        return std::string("Star(") + role_->as_str() + ")";
    }
};

class InverseRole : public Role {
  protected:
    const Role *role_;

  public:
    explicit InverseRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    ~InverseRole() override = default;
    const Role* clone() const override {
        return new InverseRole(role_);
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

    std::string as_str() const override {
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
    const Role* clone() const override {
        return new RoleRestriction(role_, restriction_);
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

    std::string as_str() const override {
        return std::string("Restriction(") + role_->as_str() + "," + restriction_->as_str() + ")";
    }
};

class Feature {
  public:
    Feature() = default;
    virtual ~Feature() = default;
    virtual const Feature* clone() const = 0;

    virtual int complexity() const = 0;
    virtual int value(const Cache &cache, const Sample &sample, const State &state) const = 0;
    virtual std::string as_str() const = 0;

    std::string as_str_with_complexity() const {
        return std::to_string(complexity()) + "." + as_str();
    }

    friend std::ostream& operator<<(std::ostream &os, const Feature &f) {
        return os << f.as_str() << std::flush;
    }

    virtual bool is_boolean() const = 0;
};

class NullaryAtomFeature : public Feature {
protected:
    const Predicate* predicate_;

public:
    explicit NullaryAtomFeature(const Predicate* predicate) : predicate_(predicate) { }
    ~NullaryAtomFeature() override = default;

    const Feature* clone() const override {
        return new NullaryAtomFeature(*this);
    }

    int complexity() const override { // Nullary atoms have complexity 0 by definition
        return 1;
    }

    int value(const Cache &cache, const Sample &sample, const State &state) const override {
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
    std::string as_str() const override {
        return std::string("Atom[") + predicate_->name_ + "]";
    }

    bool is_boolean() const override { return true; }
};

class BooleanFeature : public Feature {
  protected:
    const Concept *concept_;

  public:
    explicit BooleanFeature(const Concept *concept) : Feature(), concept_(concept) { }
    ~BooleanFeature() override = default;
    const Feature* clone() const override {
        return new BooleanFeature(concept_);
    }

    int complexity() const override {
        return concept_->complexity();
    }
    int value(const Cache &cache, const Sample &sample, const State &state) const override {
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
    std::string as_str() const override {
        return std::string("Boolean[") + concept_->as_str() + "]";
    }

    bool is_boolean() const override { return true; }
};

class NumericalFeature : public Feature {
  protected:
    const Concept *concept_;

  public:
    explicit NumericalFeature(const Concept *concept) : Feature(), concept_(concept) { }
    ~NumericalFeature() override = default;
    const Feature* clone() const override {
        return new NumericalFeature(concept_);
    }

    int complexity() const override {
        return concept_->complexity();
    }
    int value(const Cache &cache, const Sample &sample, const State &state) const override {
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
    std::string as_str() const override {
        return std::string("Numerical[") + concept_->as_str() + "]";
    }

    bool is_boolean() const override { return false; }
};

class DistanceFeature : public Feature {
  protected:
    const Concept *start_;
    const Concept *end_;
    const RoleRestriction *role_;

    mutable bool valid_cache_;
    mutable std::vector<int> cached_distances_;
    mutable bool denotation_is_constant_;
    mutable bool all_values_greater_than_zero_;

  public:
    DistanceFeature(const Concept *start, const Concept *end, const RoleRestriction *role)
      : Feature(),
        start_(start),
        end_(end),
        role_(role),
        valid_cache_(false),
        denotation_is_constant_(false),
        all_values_greater_than_zero_(false) {
    }
    ~DistanceFeature() override = default;
    const Feature* clone() const override {
        DistanceFeature *f = new DistanceFeature(start_, end_, role_);
        f->valid_cache_ = valid_cache_;
        f->cached_distances_ = cached_distances_;
        f->denotation_is_constant_ = denotation_is_constant_;
        f->all_values_greater_than_zero_ = all_values_greater_than_zero_;
        return f;
    }

    int complexity() const override {
        return 1 + start_->complexity() + end_->complexity() + role_->complexity();
    }
    int value(const Cache &cache, const Sample &sample, const State &state) const override {
        assert(sample.state(state.id()).id() == state.id());
        if( !valid_cache_ ) {
            const sample_denotation_t *start_d = cache.find_sample_denotation(start_->as_str());
            assert((start_d != nullptr) && (start_d->size() == sample.num_states()));
            const sample_denotation_t *end_d = cache.find_sample_denotation(end_->as_str());
            assert((end_d != nullptr) && (end_d->size() == sample.num_states()));
            const sample_denotation_t *role_d = cache.find_sample_denotation(role_->as_str());
            assert((role_d != nullptr) && (role_d->size() == sample.num_states()));

            denotation_is_constant_ = true;
            all_values_greater_than_zero_ = true;
            cached_distances_ = std::vector<int>(sample.num_states(), std::numeric_limits<int>::max());
            for( int i = 0; i < sample.num_states(); ++i ) {
                const State &state = sample.state(i);
                assert(state.id() == i);
                const state_denotation_t *start_sd = (*start_d)[i];
                assert((start_sd != nullptr) && (start_sd->size() == state.num_objects()));
                const state_denotation_t *end_sd = (*end_d)[i];
                assert((end_sd != nullptr) && (end_sd->size() == state.num_objects()));
                const state_denotation_t *role_sd = (*role_d)[i];
                assert((role_sd != nullptr) && (role_sd->size() == state.num_objects() * state.num_objects()));
                int distance = compute_distance(state.num_objects(), *start_sd, *end_sd, *role_sd);
                denotation_is_constant_ = denotation_is_constant_ && ((i == 0) || (cached_distances_[i - 1] == distance));
                all_values_greater_than_zero_ = all_values_greater_than_zero_ && (distance > 0);
                cached_distances_[i] = distance;
            }
            valid_cache_ = true;
        }
        return cached_distances_[state.id()];
    }
    std::string as_str() const override {
        return std::string("Distance[") + start_->as_str() + ";" + role_->as_str() + ";" + end_->as_str() + "]";
    }

    bool valid_cache() const {
        return valid_cache_;
    }
    bool denotation_is_constant() const {
        return denotation_is_constant_;
    }
    bool all_values_greater_than_zero() const {
        return all_values_greater_than_zero_;
    }

    int compute_distance(int num_objects,
                         const state_denotation_t &start_sd,
                         const state_denotation_t &end_sd,
                         const state_denotation_t &role_sd) const {
        // create adjacency lists
        std::vector<std::vector<int> > adj(num_objects);
        for( int i = 0; i < num_objects; ++i ) {
            for( int j = 0; j < num_objects; ++j ) {
                if( role_sd[i * num_objects + j] )
                    adj[i].emplace_back(j);
            }
        }

        // locate start vertex
        int start = -1;
        for( int i = 0; i < num_objects; ++i ) {
            if( start_sd[i] ) {
                start = i;
                break;
            }
        }
        assert(start != -1);

        // check whether distance is 0
        if( end_sd[start] ) return 0;

        // apply breadth-first search from start vertex
        std::vector<int> distances(num_objects, -1);
        distances[start] = 0;

        std::deque<std::pair<int, int> > q;
        for( int i = 0; i < adj[start].size(); ++i ) {
            if( end_sd[adj[start][i]] )
                return 1;
            else
                q.emplace_back(adj[start][i], 1);
        }

        while( !q.empty() ) {
            std::pair<int, int> p = q.front();
            q.pop_front();
            if( distances[p.first] == -1 ) {
                distances[p.first] = p.second;
                for( int i = 0; i < adj[p.first].size(); ++i ) {
                    if( end_sd[adj[p.first][i]] )
                        return 1 + p.second;
                    else
                        q.emplace_back(adj[p.first][i], 1 + p.second);
                }
            }
        }
        return std::numeric_limits<int>::max();
    }

    bool is_boolean() const override { return false; }
};

class Factory {
  protected:
    const std::string name_;
    const std::vector<std::string> goal_concepts_;
    const std::vector<std::string> nominals_;
    std::vector<const Role*> basis_roles_;
    std::vector<const Concept*> basis_concepts_;

    int complexity_bound_;

    mutable std::vector<const Role*> roles_;

    // A layered set of concepts, concepts_[k] contains all concepts generated in the k-th application of
    // the concept grammar
    mutable std::vector<std::vector<const Concept*> > concepts_;
    mutable std::vector<const Feature*> features_;

  public:
    Factory(std::string name, const std::vector<std::string>& nominals, int complexity_bound)
      : name_(std::move(name)),
//        goal_concepts_(goal_concepts),
        nominals_(nominals),
        complexity_bound_(complexity_bound) {
    }
    virtual ~Factory() = default;

    const std::string& name() const {
        return name_;
    }

    void insert_basis(const Role *role) {
        basis_roles_.push_back(role);
    }
    void insert_basis(const Concept *concept) {
        basis_concepts_.push_back(concept);
    }
    void set_complexity_bound(int complexity_bound) {
        complexity_bound_ = complexity_bound;
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

    void insert_new_denotation_by_name(Cache &cache, const std::string &name, const sample_denotation_t *d) const {
        cache.find_or_insert_sample_denotation_by_name(name, *d);
    }
    int insert_new_concept(Cache &cache, const Concept *concept, const sample_denotation_t *d) const {
        assert(!concepts_.empty());
        concepts_.back().push_back(concept);
        cache.find_or_insert_sample_denotation(*d, concept->as_str());
        return concepts_.back().size();
    }
    int insert_new_role(Cache &cache, const Role *role, const sample_denotation_t *d) const {
        roles_.push_back(role);
        cache.find_or_insert_sample_denotation(*d, role->as_str());
        return roles_.size();
    }

    bool is_superfluous(Cache &cache, const sample_denotation_t *d) const {
        return cache.find_sample_denotation(*d) != nullptr;
    }
    std::pair<bool, const sample_denotation_t*>
      is_superfluous_or_exceeds_complexity_bound(const Base &base, Cache &cache, const Sample &sample, bool do_not_check_complexity = false) const {
        if( !do_not_check_complexity && (base.complexity() > complexity_bound_) ) {
            return std::make_pair(true, nullptr);
        } else {
            const sample_denotation_t *d = base.denotation(cache, sample, false);
            return std::make_pair(is_superfluous(cache, d), d);
        }
    }

    // apply one iteration of the concept generation grammar
    // new concepts are left on a new last layer in concepts_
    // new concepts are non-redundant if sample != nullptr
    int advance_step(Cache &cache, const Sample &sample) const {
        int num_pruned_concepts = 0;
        if( concepts_.empty() ) {
            concepts_.emplace_back();
            for (auto concept : basis_concepts_) {
                std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(*concept, cache, sample);
                if( !p.first ) insert_new_concept(cache, concept->clone(), p.second);
                //else std::cout << "PRUNE: " + concept.as_str() << std::endl;
                delete p.second;
                num_pruned_concepts += p.first;
            }

        } else {
            bool is_first_non_basis_iteration = (concepts_.size() == 1);
            // extract concepts in existing concept layers
            std::vector<const Concept*> last_concept_layer = concepts_.back();
            std::vector<const Concept*> concepts_in_all_layers_except_last;
            for( int layer = 0; layer < int(concepts_.size()); ++layer )
                concepts_in_all_layers_except_last.insert(concepts_in_all_layers_except_last.end(), concepts_[layer].begin(), concepts_[layer].end());

            // create new concept layer
            concepts_.emplace_back();


            // Insert the negation of all primitive concepts only (not applied in later iterations)
            if (is_first_non_basis_iteration) {
                for (auto concept:last_concept_layer) {
                    NotConcept not_concept(concept);
                    auto p = is_superfluous_or_exceeds_complexity_bound(not_concept, cache, sample);
                    if (!p.first) insert_new_concept(cache, not_concept.clone(), p.second);
                    //else std::cout << "PRUNE: " + not_concept.as_str() << std::endl;
                    delete p.second;
                    num_pruned_concepts += p.first;
                }
            }

            // generate exist, and forall for concepts in last layer
            for( int i = 0; i < int(last_concept_layer.size()); ++i ) {
                const Concept *concept = last_concept_layer[i];
#if 0 // Don't negate in layers > 0
                NotConcept not_concept(concept);
                std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(not_concept, cache, sample);
                if( !p.first ) insert_new_concept(cache, not_concept.clone(), p.second);
                //else std::cout << "PRUNE: " + not_concept.as_str() << std::endl;
                delete p.second;
                num_pruned_concepts += p.first;
#endif

                // TODO: Merge the two loops below (for exist- and forall- concepts) into one single loop.
                // TODO: We are keeping the separate ATM to be able to compare better with the Python version
                for (auto role : roles_) {
                    ExistsConcept exists_concept(concept, role);
                    std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(exists_concept, cache, sample);
                    if( !p.first ) insert_new_concept(cache, exists_concept.clone(), p.second);
                    //else std::cout << "PRUNE: " + exists_concept.as_str() << std::endl;
                    delete p.second;
                    num_pruned_concepts += p.first;
                }

                for (auto role : roles_) {
                    ForallConcept forall_concept(concept, role);
                    std::pair<bool, const sample_denotation_t*> q = is_superfluous_or_exceeds_complexity_bound(forall_concept, cache, sample);
                    if( !q.first ) insert_new_concept(cache, forall_concept.clone(), q.second);
                    //else std::cout << "PRUNE: " + forall_concept.as_str() << std::endl;
                    delete q.second;
                    num_pruned_concepts += q.first;
                }
            }

            // Insert equal concepts based on the already-fixed set of roles. We only do this once, and it'd thus be
            // more sensible to put this code elsewhere, but we want to generate concepts in the same order than
            // the python code for easier comparability
            if (is_first_non_basis_iteration) {
                for (unsigned i = 0; i < roles_.size(); ++i) {
                    for (unsigned j = i + 1; j < roles_.size(); ++j) {
                        EqualConcept eq_concept(roles_[i], roles_[j]);
                        auto p = is_superfluous_or_exceeds_complexity_bound(eq_concept, cache, sample);
                        if (!p.first) insert_new_concept(cache, eq_concept.clone(), p.second);
                        delete p.second;
                        num_pruned_concepts += p.first;
                    }
                }
            }

            // generate conjunction of concepts in last layer with all other
            // concepts, avoiding redudant pairs
            for( int i = 0; i < int(last_concept_layer.size()); ++i ) {
                const Concept *concept1 = last_concept_layer[i];

                for( int j = 0; j < int(concepts_in_all_layers_except_last.size()); ++j ) {
                    const Concept *concept2 = concepts_in_all_layers_except_last[j];
                    AndConcept and_concept(concept1, concept2);
                    std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(and_concept, cache, sample);
                    if( !p.first ) insert_new_concept(cache, and_concept.clone(), p.second);
                    //else std::cout << "PRUNE: " + and_concept.as_str() << std::endl;
                    delete p.second;
                    num_pruned_concepts += p.first;
                }

                for( int j = 1 + i; j < int(last_concept_layer.size()); ++j ) {
                    const Concept *concept2 = last_concept_layer[j];
                    AndConcept and_concept(concept1, concept2);
                    std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(and_concept, cache, sample);
                    if( !p.first ) insert_new_concept(cache, and_concept.clone(), p.second);
                    //else std::cout << "PRUNE: " + and_concept.as_str() << std::endl;
                    delete p.second;
                    num_pruned_concepts += p.first;
                }
            }
        }
        return num_pruned_concepts;
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
        std::cout << "BASIS: #concepts=" << basis_concepts_.size()
                  << ", #roles=" << basis_roles_.size()
                  << std::endl;
    }

    void generate_goal_concepts(const Sample &sample) {
        /*
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

        for (const auto& gc:goal_concepts_) {
            if (sample.predicates())
            insert_basis(new NominalConcept(nominal));
        }
        */
    }

    int generate_roles(Cache &cache, const Sample &sample) const {
        assert(roles_.empty());
        for( int i = 0; i < int(basis_roles_.size()); ++i ) {
            const Role &role = *basis_roles_[i];
            std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(role, cache, sample);
            if( !p.first ) {
                insert_new_role(cache, role.clone(), p.second);
            } else {
                // we cannot just obviate this role since another
                // role denotation may depend on this denotation
                //std::cout << "PRUNE: " + role.as_str() << std::endl;
                insert_new_denotation_by_name(cache, role.as_str(), p.second);
            }
            delete p.second;
        }

        for( int i = 0; i < int(basis_roles_.size()); ++i ) {
            const Role *role = basis_roles_[i];

            PlusRole plus_role(role);
            std::pair<bool, const sample_denotation_t*> p1 = is_superfluous_or_exceeds_complexity_bound(plus_role, cache, sample);
            if( !p1.first ) {
                insert_new_role(cache, plus_role.clone(), p1.second);
            } else {
                // we cannot just obviate this role since another
                // role denotation may depend on this denotation
                //std::cout << "PRUNE: " + plus_role.as_str() << std::endl;
                insert_new_denotation_by_name(cache, plus_role.as_str(), p1.second);
            }
            delete p1.second;
#if 0  // ATM we deactivate these roles
            StarRole star_role(role);
            std::pair<bool, const sample_denotation_t*> p2 = is_superfluous_or_exceeds_complexity_bound(star_role, cache, sample);
            if( !p2.first ) {
                insert_new_role(cache, star_role.clone(), p2.second);
            } else {
                // we cannot just obviate this role since another
                // role denotation may depend on this denotation
                //std::cout << "PRUNE: " + star_role.as_str() << std::endl;
                insert_new_denotation_by_name(cache, star_role.as_str(), p2.second);
            }
            delete p2.second;
#endif
            InverseRole inverse_role(role);
            std::pair<bool, const sample_denotation_t*> p3 = is_superfluous_or_exceeds_complexity_bound(inverse_role, cache, sample);
            if( !p3.first ) {
                const Role *clone = inverse_role.clone();
                insert_new_role(cache, clone, p3.second);

                PlusRole plus_inverse_role(clone);
                std::pair<bool, const sample_denotation_t*> p = is_superfluous_or_exceeds_complexity_bound(plus_inverse_role, cache, sample);
                if( !p.first ) {
                    insert_new_role(cache, plus_inverse_role.clone(), p.second);
                } else {
                    // we cannot just obviate this role since another
                    // role denotation may depend on this denotation
                    //std::cout << "PRUNE: " + plus_inverse_role.as_str() << std::endl;
                    insert_new_denotation_by_name(cache, plus_inverse_role.as_str(), p.second);
                }
                delete p.second;

#if 0  // ATM we deactivate these roles

                PlusRole star_inverse_role(clone);
                std::pair<bool, const sample_denotation_t*> q = is_superfluous_or_exceeds_complexity_bound(star_inverse_role, cache, sample);
                if( !q.first ) {
                    insert_new_role(cache, star_inverse_role.clone(), q.second);
                } else {
                    // we cannot just obviate this role since another
                    // role denotation may depend on this denotation
                    //std::cout << "PRUNE: " + star_inverse_role.as_str() << std::endl;
                    insert_new_denotation_by_name(cache, star_inverse_role.as_str(), q.second);
                }
                delete q.second;
#endif

            } else {
                // we cannot just obviate this role since another
                // role denotation may depend on this denotation
                //std::cout << "PRUNE: " + inverse_role.as_str() << std::endl;
                insert_new_denotation_by_name(cache, inverse_role.as_str(), p3.second);
            }
            delete p3.second;
        }
        std::cout << "ROLES: #roles=" << roles_.size() << std::endl;
        return roles_.size();
    }

    std::vector<const Concept*> generate_concepts(Cache &cache, const Sample &sample) const {
        int num_concepts = 0;
        bool some_new_concepts = true;
        for( int iteration = 0; some_new_concepts; ++iteration ) {
            std::cout << "DL::Factory: iteration=" << iteration
                      << ", #concepts=" << num_concepts
                      << ", #concepts-in-last-layer=" << (concepts_.empty() ? 0 : concepts_.back().size())
                      << std::endl;

            int num_concepts_before_step = num_concepts;
            int num_pruned_concepts = advance_step(cache, sample);
            num_concepts += concepts_.empty() ? 0 : concepts_.back().size();
            some_new_concepts = num_concepts > num_concepts_before_step;

            std::cout << "DL::Factory: advance-step:"
                      << " #concepts-in-layer=" << concepts_.back().size()
                      << ", #pruned-concepts=" << num_pruned_concepts
                      << std::endl;
        }
        assert(concepts_.back().empty());
        concepts_.pop_back();

        // Important to use this vector from here on, as it is sorted by complexity
        auto all_concepts = consolidate_concepts();
        assert(all_concepts.size() == num_concepts);

        std::cout << "DL::Factory: name=" << name_ << ", #concepts-final=" << num_concepts << std::endl;
        return all_concepts;
    }

    int generate_features(const std::vector<const Concept*>& concepts, Cache &cache, const Sample &sample) const {
        using cache_t = std::unordered_map<feature_sample_denotation_t, const Feature*, utils::container_hash<feature_sample_denotation_t> >;
        cache_t seen_denotations;
        unsigned num_nullary_features = 0, num_boolean_features = 0, num_numeric_features = 0, num_distance_features = 0;

        // create features that derive from nullary predicates
        for (const auto& predicate:sample.predicates()) {
            if (predicate.arity() == 0) {
                features_.push_back(new NullaryAtomFeature(&predicate));
                num_nullary_features++;
            }
        }

        // create boolean/numerical features from concepts
        for (const Concept* c:concepts) {
            //std::cout << c->as_str() << std::endl;
            const sample_denotation_t *d = cache.find_sample_denotation(c->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));

            // generate feature denotation associated to concept's sample denotation
            feature_sample_denotation_t fd;
            fd.reserve(sample.num_states());

            // track whether feature would be boolean or numeric, or non-informative
            bool boolean_feature = true;
            bool denotation_is_constant = true;
            bool all_values_greater_than_zero = true;
            int previous_value = -1;

            for( int j = 0; j < sample.num_states(); ++j ) {
                const State &state = sample.state(j);
                assert(state.id() == j);
                const state_denotation_t *sd = (*d)[j];
                assert((sd != nullptr) && (sd->size() == state.num_objects()));
                int value = sd->cardinality();
                boolean_feature = boolean_feature && (value < 2);
                denotation_is_constant = (previous_value == -1) || (denotation_is_constant && (previous_value == value));
                all_values_greater_than_zero = all_values_greater_than_zero && (value > 0);
                previous_value = value;
                fd.push_back(value);
            }

            if( !denotation_is_constant && !all_values_greater_than_zero ) {
                cache_t::const_iterator it = seen_denotations.find(fd);
                if( it == seen_denotations.end() ) { // The feature denotation is new, keep the feature
                    num_boolean_features += boolean_feature;
                    num_numeric_features += !boolean_feature;
                    const Feature *feature = boolean_feature ? static_cast<Feature*>(new BooleanFeature(c)) :
                            static_cast<Feature*>(new NumericalFeature(c));
                    features_.emplace_back(feature);
                    seen_denotations.emplace(fd, feature);
                    //std::cout << "ACCEPT: " << feature->as_str_with_complexity() << std::endl;
                } else {
                    // Make sure that we are not pruning a feature of lower complexity in favor of a feature of higher
                    // complexity!
                    assert(it->second->complexity() <= c->complexity());
//                    std::cout << "REJECT: " << c->as_str() << std::endl;
//                    std::cout << "PRUNED-BY: " << it->second->as_str() << std::endl;
                }
            }
        }

        // create distance features - temporarily deactivated
         num_distance_features = generate_distance_features(concepts, cache, sample);

        assert(num_nullary_features+num_boolean_features+num_numeric_features+num_distance_features == features_.size());
        std::cout << "FEATURES: #features=" << features_.size()
                  << ", #nullary="   << num_nullary_features
                  << ", #boolean="   << num_boolean_features
                  << ", #numerical=" << num_numeric_features
                  << ", #distance="  << num_distance_features
                  << std::endl;

        return features_.size();
    }

    int generate_distance_features(const std::vector<const Concept*>& concepts, Cache &cache, const Sample &sample) const {
        using cache_t = std::unordered_map<feature_sample_denotation_t, const Feature*, utils::container_hash<feature_sample_denotation_t> >;
        cache_t seen_denotations;

        // identify concepts with singleton denotation across all states
        // as these are the candidates for start concepts
        std::vector<const Concept*> start_concepts;
        for (const Concept* c:concepts) {
            const sample_denotation_t *d = cache.find_sample_denotation(c->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));
            bool singleton_denotations = true;
            for( int j = 0; singleton_denotations && (j < sample.num_states()); ++j )
                singleton_denotations = (*d)[j]->cardinality() == 1;
            if( singleton_denotations ) {
                //std::cout << "GOOD (START) CONCEPT: " << c.as_str_with_complexity() << std::endl;
                start_concepts.push_back(c);
            }
        }

        // create role restrictions to be used in distance features
        Cache cache_for_role_restrictions;
        std::vector<const RoleRestriction*> role_restrictions;
        for (const Role* r:roles_) {
            for (const Concept* c:concepts) {
                RoleRestriction role_restriction(r, c);
                if( 3 + role_restriction.complexity() > complexity_bound_ ) continue;
                const sample_denotation_t *d = role_restriction.denotation(cache, sample, false);
                if( !is_superfluous(cache_for_role_restrictions, d) ) {
                    assert(cache_for_role_restrictions.cache1().find(d) == cache_for_role_restrictions.cache1().end());
                    role_restrictions.push_back(static_cast<const RoleRestriction*>(role_restriction.clone()));
                    cache_for_role_restrictions.find_or_insert_sample_denotation(*d, role_restriction.as_str());
                    assert(cache_for_role_restrictions.cache1().find(d) != cache_for_role_restrictions.cache1().end());
                    //std::cout << "ACCEPT RR(sz=" << cache_for_role_restrictions.cache1().size() << "): " + role_restriction.as_str_with_complexity() << std::endl;
                } else {
                    //std::cout << "PRUNE RR: " + role_restriction.as_str() << std::endl;
                }
                delete d;
            }
        }

        // create distance features
        int num_distance_features = 0;
        for( int i = 0; i < int(start_concepts.size()); ++i ) {
            const Concept &start = *start_concepts[i];
            for (const Concept* end:concepts) {
                for( int k = 0; k < int(role_restrictions.size()); ++k ) {
                    const RoleRestriction &role = *role_restrictions[k];
                    DistanceFeature df(&start, end, &role);
                    if( df.complexity() > complexity_bound_ ) continue;
                    //std::cout << "TESTING: " << df.as_str_with_complexity() << std::endl;

                    // fill cache with denotations for start and end concepts
                    const sample_denotation_t *ds = cache.find_sample_denotation(start.as_str());
                    assert((ds != nullptr) && (ds->size() == sample.num_states()));
                    cache_for_role_restrictions.find_or_insert_sample_denotation(*ds, start.as_str());
                    const sample_denotation_t *de = cache.find_sample_denotation(end->as_str());
                    assert((de != nullptr) && (de->size() == sample.num_states()));
                    cache_for_role_restrictions.find_or_insert_sample_denotation(*de, end->as_str());

                    // check whether df is novel
                    feature_sample_denotation_t fd(sample.num_states());
                    for( int index = 0; index < sample.num_states(); ++index ) {
                        assert(sample.state(index).id() == index);
                        const State &state = sample.state(index);
                        int value = df.value(cache_for_role_restrictions, sample, state);
                        fd[index] = value;
                    }

                    cache_t::const_iterator it = seen_denotations.find(fd);
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
        return num_distance_features;;
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

        os << "All Roles under complexity " << complexity_bound_ << "(sz=" << roles_.size() << "): ";
        for( int i = 0; i < int(roles_.size()); ++i ) {
            os << roles_[i]->as_str_with_complexity();
            if( 1 + i < int(roles_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "All concepts (by layer) under complexity " << complexity_bound_ << ": " << std::endl;
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
        for (auto &feature : features_) {
            if( dynamic_cast<const NullaryAtomFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Boolean features: ";
        need_comma = false;
        for (auto &feature : features_) {
            if( dynamic_cast<const BooleanFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Numerical features: ";
        need_comma = false;
        for (auto &feature : features_) {
            if( dynamic_cast<const NumericalFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Distance features: ";
        need_comma = false;
        for (auto &feature : features_) {
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
        int num_features = features_.size();
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
        int num_features = features_.size();

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
    }

    void log_all_concepts_and_features(const std::vector<const Concept*>& concepts,
            const SLTP::DL::Cache &cache, const SLTP::DL::Sample &sample,
            const std::string& workspace);

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

}; // DL namespace

}; // SLTP namespace

#endif

