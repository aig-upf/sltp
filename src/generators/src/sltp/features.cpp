
#include <sltp/features.hxx>
#include <blai/transitions.h>

#include <fstream>
#include <common/helpers.h>


using namespace std;

namespace sltp::dl {


const state_denotation_t& Cache::retrieveDLDenotation(
        const DLBaseElement& element, const State &state, std::size_t expected_size) const {
    const sample_denotation_t& sd = find_sample_denotation(element, expected_size);
    return *sd[state.id()];
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

        const auto m = sample.num_states();
        for (unsigned i = 0; i < m; ++i) {
            const State &state = sample.state(i);
            const auto& oidx = sample.instance(i).object_index();

            for (const Concept *c:concepts) {
                const state_denotation_t &denotation = cache.retrieveDLDenotation(*c, state, m);
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

        for (unsigned i = 0; i < m; ++i) {
            const State &state = sample.state(i);
            const auto& oidx = sample.instance(i).object_index();
            unsigned n = sample.num_objects(i);


            for (const Role *r:roles_) {
                const state_denotation_t &denotation = cache.retrieveDLDenotation(*r, state, m);
                of << "s_" << i << "[" << r->as_str() << "] = {";
                bool need_comma = false;

                for (unsigned idx = 0; idx < denotation.size(); ++idx) {
                    if (denotation[idx]) {
                        if (need_comma) of << ", ";
                        unsigned o1 = idx / n;
                        unsigned o2 = idx % n;
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

        for (unsigned i = 0; i < m; ++i) {
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
        if (attempt_cardinality_feature_insertion(c, cache, sample, transitions, seen_denotations, false)) {
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
        attempt_cardinality_feature_insertion(c, cache, sample, transitions, seen_denotations, true);
    }

    // create comparison features here so that only cardinality features are used to build them
    generate_comparison_features(features_, cache, sample, transitions, seen_denotations);

    // create distance features
    generate_distance_features(concepts, cache, sample, transitions, seen_denotations);

    // create conditional features from boolean conditions and numeric bodies
    generate_conditional_features(features_, cache, sample, transitions, seen_denotations);

    print_feature_count();
}

bool Factory::attempt_cardinality_feature_insertion(
        const Concept* c,
        Cache &cache,
        const Sample &sample,
        const TransitionSample& transitions,
        feature_cache_t& seen_denotations,
        bool check_redundancy)
{
    // TODO We could compute the whole denotation with one single fetch of the underlying concept
    //      denotation. By doing as below, `compute_feature_sample_denotation` calls n times
    //      that same fetch, one per state in the sample.
    NumericalFeature nf(c);
    SampleDenotationProperties properties;
    auto fd = compute_feature_sample_denotation(nf, sample, cache, properties);

    if (prune_feature_denotation(nf, fd, properties, sample, transitions, seen_denotations, check_redundancy)) {
        return false;
    }

    const Feature *feature = properties.denotation_is_bool ?
            static_cast<Feature*>(new BooleanFeature(c)) : static_cast<Feature*>(new NumericalFeature(c));

    features_.emplace_back(feature);
    seen_denotations.emplace(fd, feature);
    return true;
}

bool Factory::check_some_transition_pair_distinguished(const feature_sample_denotation_t &fsd, const Sample &sample,
                                                       const TransitionSample &transitions) {
// Make sure that the feature is useful to distinguish at least some pair of transitions
// coming from the same instance.
// Since the notion of distinguishability is transitive, we'll just check for one pair distinguished from the previous
// one
        int prev_instance_id = -1;
        bool feature_can_distinguish_some_transition = false;
        int last_sf = -1, last_sfprime = -1;
        for (auto s:transitions.all_alive()) {
            const State &state = sample.state(s);

            int sf = fsd[s];

            for (unsigned sprime:transitions.successors(s)) {
                int sfprime = fsd[sprime];

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
        return feature_can_distinguish_some_transition;
    }

// TODO Ideally we'd want to unify this function with Factory::attempt_cardinality_feature_insertion, but
//  at we're not there yet
bool Factory::attempt_feature_insertion(
        const Feature* feature,
        unsigned bound,
        Cache &cache,
        const Sample &sample,
        const TransitionSample& transitions,
        feature_cache_t& seen_denotations,
        bool check_redundancy)
{
    if (feature->complexity() > bound) return false;

    SampleDenotationProperties properties;
    auto fd = compute_feature_sample_denotation(*feature, sample, cache, properties);

    if (prune_feature_denotation(
            *feature, fd, properties, sample, transitions, seen_denotations, check_redundancy)) {
        return false;
    }

    features_.emplace_back(feature);
    seen_denotations.emplace(fd, feature);

    return true;
}

bool Factory::prune_feature_denotation(
        const Feature& f,
        const feature_sample_denotation_t& fd,
        const SampleDenotationProperties& properties,
        const Sample &sample,
        const TransitionSample& transitions,
        feature_cache_t& seen_denotations,
        bool check_redundancy)
{
    // We want to determine:
    // - whether the feature is boolean or numeric (this is simply determined empirically: we consider it boolean
    //   if its value is always 0 or 1, and numeric otherwise),
    // - whether the full sample denotation of the feature coincides with some previous feature and hence we can prune
    //   it,
    // - whether the feature is not truly informative for our encodings, e.g. because it has the same variation over all
    //   transitions in the sample, or similar
    if (!check_redundancy) return false;

    if (properties.denotation_is_constant) return true;

    auto it = seen_denotations.find(fd);
    bool is_new = (it == seen_denotations.end());

    if (!is_new) {
//         std::cout << "REJECT (Redundant): " << f.as_str_with_complexity() << std::endl;
        // Make sure that we don't prune a feature of lower complexity in favor of a feature of higher complexity
        // This should come for free, since features are ordered in increasing complexity
        if (it->second->complexity() > f.complexity()) {
            std::cout << Utils::warning()
                      <<  "Feature " + f.as_str_with_complexity() + " has been pruned in favor of more complex "
                          + it->second->as_str_with_complexity() << std::endl;
        }

        return true;
    }

    if (!check_some_transition_pair_distinguished(fd, sample, transitions)) {
//         std::cout << "REJECT (NO DISTINCTION): " << f.as_str_with_complexity() << std::endl;
        return true;
    }

//    std::cout << "ACCEPT: " << f.as_str_with_complexity() << std::endl;
    return false;
}

feature_sample_denotation_t Factory::compute_feature_sample_denotation(
        const Feature& feature, const Sample &sample, const Cache &cache, SampleDenotationProperties& properties) {

    const auto m = sample.num_states();

    feature_sample_denotation_t fd;
    fd.reserve(m);

    properties.denotation_is_bool = true;
    properties.denotation_is_constant = true;
    int previous_value = -1;

    for (unsigned sid = 0; sid < m; ++sid) {
        const State &state = sample.state(sid);
        assert(state.id() == sid);

        int value = feature.value(cache, sample, state);
        fd.push_back(value);

        properties.denotation_is_bool = properties.denotation_is_bool && (value < 2);
        properties.denotation_is_constant = (previous_value == -1)
                                            || (properties.denotation_is_constant && (previous_value == value));
        previous_value = value;
    }

    return fd;
}

const sample_denotation_t&
Cache::find_sample_denotation(const DLBaseElement& element, std::size_t expected_size) const {
    auto it = cache2_.find(element.id());
    assert (it != cache2_.end());
    assert(it->second != nullptr);
    assert(it->second->size() == expected_size);
    return *it->second;
}

unsigned long DLBaseElement::global_id = 0;

} // namespaces
