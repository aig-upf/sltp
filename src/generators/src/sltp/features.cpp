
#include <sltp/features.hxx>
#include <blai/transitions.h>

#include <fstream>
#include <common/helpers.h>


using namespace std;

namespace sltp::dl {


const state_denotation_t& Cache::retrieve_concept_denotation(const Concept &element, const State &state) const {
    const sample_denotation_t *d = find_sample_denotation(element.as_str());
    assert(d);
    const state_denotation_t *sd = (*d)[state.id()];
    return *sd;
}


const state_denotation_t& Cache::retrieve_role_denotation(const Role &element, const State &state) const {
    const sample_denotation_t *d = find_sample_denotation(element.as_str());
    assert(d);
    const state_denotation_t *sd = (*d)[state.id()];
    return *sd;
}

void Factory::log_all_concepts_and_features(
        const std::vector<const Concept*>& concepts,
        const Cache &cache, const Sample &sample, const string &workspace,
        bool print_denotations) {

    if (print_denotations) {
        std::cout << "Printing concept, role and feature denotations to " << workspace
                  << "/*-denotations.io.txt" << std::endl;
        // Print concept denotations
        std::string output(workspace + "/concept-denotations.io.txt");
        std::ofstream of(output);
        if (of.fail()) throw std::runtime_error("Could not open filename '" + output + "'");

        for (unsigned i = 0; i < sample.num_states(); ++i) {
            const State &state = sample.state(i);
            const auto& oidx = sample.instance(i).object_index();

            for (const Concept *c:concepts) {
                const state_denotation_t &denotation = cache.retrieve_concept_denotation(*c, state);
                of << "s_" << i << "[" << c->as_str() << "] = {";
                bool need_comma = false;
                for (unsigned atom = 0; atom < denotation.size(); ++atom) {
                    if (denotation[atom]) {
                        if (need_comma) of << ", ";
                        of << oidx.right.at(atom);
                        need_comma = true;
                    }
                }
                of << "}" << std::endl;
            }
        }
        of.close();


        // Print role denotations
        output = workspace + "/role-denotations.io.txt";
        of = std::ofstream(output);
        if (of.fail()) throw std::runtime_error("Could not open filename '" + output + "'");

        for (unsigned i = 0; i < sample.num_states(); ++i) {
            const State &state = sample.state(i);
            const auto& oidx = sample.instance(i).object_index();
            unsigned m = sample.num_objects(i);


            for (const Role *r:roles_) {
                const state_denotation_t &denotation = cache.retrieve_role_denotation(*r, state);
                of << "s_" << i << "[" << r->as_str() << "] = {";
                bool need_comma = false;

                for (unsigned idx = 0; idx < denotation.size(); ++idx) {
                    if (denotation[idx]) {
                        if (need_comma) of << ", ";
                        unsigned o1 = idx / m;
                        unsigned o2 = idx % m;
                        of << "(" << oidx.right.at(o1) << ", " << oidx.right.at(o2) << ")";
                        need_comma = true;
                    }
                }
                of << "}" << std::endl;
            }
        }
        of.close();

        // Print feature denotations
        output = workspace + "/feature-denotations.io.txt";
        of = std::ofstream(output);
        if (of.fail()) throw std::runtime_error("Could not open filename '" + output + "'");

        for (unsigned i = 0; i < sample.num_states(); ++i) {
            const State &state = sample.state(i);

            for (const Feature *f:features_) {
                of << "s_" << i << "[" << f->as_str() << "] = " << f->value(cache, sample, state) << std::endl;
            }
        }
        of.close();
    }
    // Print all generated concepts
    std::string fname1 = workspace + "/serialized-concepts.io";
    std::ofstream of1 = std::ofstream(fname1);
    if( of1.fail() ) throw std::runtime_error("Could not open filename '" + fname1 + "'");

    // Print all generated features to be unserialized from the Python frontend
    std::string fname2 = workspace + "/serialized-features.io";
    std::ofstream of2 = std::ofstream(fname2);
    if( of2.fail() ) throw std::runtime_error("Could not open filename '" + fname2 + "'");

    std::cout << "Serializing all concepts and features to:\n\t" << fname1 << "\n\t" << fname2 << std::endl;
    for (const Concept* c:concepts) {
        of1 << c->as_str() << "\t" << c->complexity() << std::endl;
    }
    of1.close();

    for (const Feature* f:features_) {
        of2 << f->as_str() << "\t" << f->complexity() << std::endl;
    }
    of2.close();
}


void Factory::generate_features(
        const std::vector<const Concept*>& concepts,
        Cache &cache, const Sample &sample,
        const TransitionSample& transitions,
        const std::vector<const Concept*>& forced_goal_features)
{
    feature_cache_t seen_denotations;

    // Insert first the features that allow us to express the goal
    // goal_features will contain the indexes of those features
    goal_features_.clear();
    for (const auto *c:forced_goal_features) {
        if (generate_cardinality_feature_if_not_redundant(c, cache, sample, transitions, seen_denotations, false)) {
            goal_features_.insert(features_.size()-1);
        }
    }
    std::cout << "A total of " << goal_features_.size() << " features were marked as goal-identifying" << std::endl;

    // create features that derive from nullary predicates
    // TODO Keep track of the denotation of nullary-atom features and prune them if necessary
    for (const auto& predicate:sample.predicates()) {
        if (predicate.arity() == 0) {
            features_.push_back(new NullaryAtomFeature(&predicate));
        }
    }

    // create boolean/numerical features from concepts
    for (const Concept* c:concepts) {
        generate_cardinality_feature_if_not_redundant(c, cache, sample, transitions, seen_denotations, true);
    }

    // create comparison features here so that only cardinality features are used to build them
    generate_comparison_features(features_, cache, sample, seen_denotations);

    // create distance features
    generate_distance_features(concepts, cache, sample, seen_denotations);

    // create conditional features from boolean conditions and numeric bodies
    generate_conditional_features(features_, cache, sample, seen_denotations);

    print_feature_count();
}

bool Factory::generate_cardinality_feature_if_not_redundant(
        const Concept* c,
        Cache &cache,
        const Sample &sample,
        const TransitionSample& transitions,
        feature_cache_t& seen_denotations,
        bool can_be_pruned)
{
    const auto m = sample.num_states();
    const sample_denotation_t *d = cache.find_sample_denotation(c->as_str());
    assert((d != nullptr) && (d->size() == m));

    // generate feature denotation associated to concept's sample denotation
    feature_sample_denotation_t fd;
    fd.reserve(m);

    // We want to determine:
    // - whether the feature is boolean or numeric (this is simply determined empirically: we consider it boolean
    //   if its value is always 0 or 1, and numeric otherwise),
    // - whether the full sample denotation of the feature coincides with some previous feature and hence we can prune
    //   it,
    // - whether the feature is not truly informative for our encodings, e.g. because it has the same variation over all
    //   transitions in the sample, or similar

    bool is_bool = true;
    for (unsigned sid = 0; sid < m; ++sid) {
        const State &state = sample.state(sid);
        assert(state.id() == sid);

        const state_denotation_t *sd = (*d)[sid];
        assert((sd != nullptr) && (sd->size() == sample.num_objects(sid)));

        int value = static_cast<int>(sd->cardinality());
        fd.push_back(value);
        is_bool = is_bool && (value < 2);
    }

    // Now the heavy checks: make sure that the feature is useful to distinguish at least some pair of transitions
    // coming from the same instance.
    // Since the notion of distinguishability is transitive, we'll just check
    int prev_instance_id = -1;
    bool feature_can_distinguish_some_transition = false;
    int last_sf = -1, last_sfprime = -1;
    for (auto s:transitions.all_alive()) {
        const State &state = sample.state(s);

        int sf = static_cast<int>((*d)[s]->cardinality());

        for (unsigned sprime:transitions.successors(s)) {
            int sfprime = static_cast<int>((*d)[sprime]->cardinality());

            if (last_sfprime > 0 && are_transitions_d1d2_distinguished(last_sf, last_sfprime, sf, sfprime)) {
                feature_can_distinguish_some_transition = true;
            }

            last_sfprime = sfprime;
            last_sf = sf;
        }

        if (prev_instance_id < 0) prev_instance_id = (int) state.instance_id();
        if (state.instance_id() != prev_instance_id) {
            last_sfprime = -1;
        }
    }


    if (can_be_pruned && !feature_can_distinguish_some_transition) return false;

    auto it = seen_denotations.find(fd);
    if( it == seen_denotations.end() ) { // The feature denotation is new, keep the feature
        const Feature *feature = is_bool ? static_cast<Feature*>(new BooleanFeature(c)) :
                                 static_cast<Feature*>(new NumericalFeature(c));
        features_.emplace_back(feature);
        seen_denotations.emplace(fd, feature);
        //std::cout << "ACCEPT: " << feature->as_str_with_complexity() << std::endl;
        return true;
    }

    // Make sure that we don't prune a feature of lower complexity in favor of a feature of higher complexity
    if (it->second->complexity() > c->complexity()) {
        throw std::runtime_error(
                "Feature " + it->second->as_str_with_complexity() + " would be pruned in favor of more complex feature "
                + c->as_str_with_complexity() + "! Check the code and fix the issue");
    }

    // From here on, the we know the feature is redundant, but still it may be that we want to force it (e.g.
    // because it is goal-identifying).

//        std::cout << "REJECT: " << c->as_str() << std::endl;
//        std::cout << "PRUNED-BY: " << it->second->as_str() << std::endl;
    if (can_be_pruned) {
        return false; // If can be pruned, we don't generate anything and return false

    } else { // Otherwise, we generate it (even if it is redundant, and return true
        const Feature *feature = is_bool ? static_cast<Feature*>(new BooleanFeature(c)) :
                                 static_cast<Feature*>(new NumericalFeature(c));
        features_.emplace_back(feature);
        return true;
    }
}

// TODO Ideally we'd want to unify this function with Factory::generate_cardinality_feature_if_not_redundant, but
//  at the moment that would impose some large performance penalty, since we would go from fetching the whole sample
//  denotation once, with an expensive call to cache.find_sample_denotation(c->as_str()), to doing that
//  once per state in the sample.
bool Factory::insert_feature_if_necessary(
        const Feature* feature, unsigned bound,
        Cache &cache, const Sample &sample, feature_cache_t& seen_denotations)
{
    if (feature->complexity() > bound) return false;

    const auto m = sample.num_states();
    feature_sample_denotation_t fd;
    fd.reserve(m);

    bool denotation_is_constant = true;
    int previous_value = -1;

    for (unsigned j = 0; j < m; ++j) {
        const State &state = sample.state(j);
        int value = feature->value(cache, sample, state);
        denotation_is_constant = (previous_value == -1) || (denotation_is_constant && (previous_value == value));
        previous_value = value;
        fd.push_back(value);
    }

    if (!denotation_is_constant) {
        if (seen_denotations.find(fd)  == seen_denotations.end()) {
            // The feature denotation is new, so let's insert it
            features_.emplace_back(feature);
            seen_denotations.emplace(fd, feature);
            return true;
        }
    }
    return false;
}

} // namespaces
