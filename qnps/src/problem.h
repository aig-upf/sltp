#ifndef PROBLEM_H
#define PROBLEM_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "basic.h"
#include "action.h"
#include "feature.h"

namespace QNP {

class Problem {
  protected:
    const std::string name_;
    std::vector<const Feature*> features_;
    std::vector<const Action*> actions_;
    std::map<std::string, const Feature*> feature_map_;

    void clone(Problem *fond, const Action *action) const;

    void create_set_actions(Problem *fond,
                            int d,
                            int loop_nesting,
                            std::vector<const Feature*> &numeric_features,
                            std::vector<const Feature*> &boolean_features,
                            std::vector<const Feature*> &bit_features,
                            std::vector<const Feature*> &q_features) const;

    void translate(Problem *fond,
                   int d,
                   int loop_nesting,
                   std::vector<const Feature*> &numeric_features,
                   std::vector<const Feature*> &boolean_features,
                   std::vector<const Feature*> &bit_features,
                   std::vector<const Feature*> &q_features) const;

  public:
    Problem(const std::string &name) : name_(name) { }
    ~Problem() {
        for( size_t i = 0; i < actions_.size(); ++i )
            delete actions_[i];
        for( size_t i = 0; i < features_.size(); ++i )
            delete features_[i];
    }

    const std::string& name() const { return name_; }
    const Feature& feature(int i) const {
        return *features_[i];
    }
    const Action& action(int i) const {
        return *actions_[i];
    }
    const Feature* feature(const std::string &name) const {
        std::map<std::string, const Feature*>::const_iterator it = feature_map_.find(name);
        return it == feature_map_.end() ? nullptr : it->second;
    }
    const std::map<std::string, const Feature*>& feature_map() const {
        return feature_map_;
    }

    void add_feature(const Feature *feature) {
        assert(feature != nullptr);
        assert(feature_map_.find(feature->name()) == feature_map_.end());
        features_.push_back(feature);
        feature_map_.insert(std::make_pair(feature->name(), features_.back()));
    }
    void add_action(const Action *action) {
        assert(action != nullptr);
        actions_.push_back(action);
    }

    Problem* translate(int d, int loop_nesting = 0) const;
#if 0
    {
        std::vector<const Feature*> numeric_features;
        std::vector<const Feature*> boolean_features;
        std::vector<const Feature*> q_features;
        std::vector<const Feature*> bit_features;
        for( size_t i = 0; i < features_.size(); ++i ) {
            if( features_[i]->numeric() )
                numeric_features.push_back(features_[i]);
            else
                boolean_features.push_back(features_[i]);
        }

        Problem *fond = new Problem(std::string("FOND_") + name_ + "_" + std::to_string(d));
        for( size_t i = 0; i < numeric_features.size(); ++i )
            fond->add_feature(new Feature(numeric_features[i]->name(), true));
        for( size_t i = 0; i < boolean_features.size(); ++i )
            fond->add_feature(new Feature(boolean_features[i]->name(), false));
        for( size_t i = 0; i < numeric_features.size(); ++i ) {
            //const Feature *feature = new Feature(std::string("Q(") + numeric_features[i]->name() + ")", false);
            const Feature *feature = new Feature(PDDL_q(numeric_features[i]->name()), false, true);
            q_features.push_back(feature);
            fond->add_feature(feature);
        }
        for( int i = 0; i <= d; ++i ) {
            const Feature *feature = new Feature(std::string("bit(") + std::to_string(i) + ")", false);
            bit_features.push_back(feature);
            fond->add_feature(feature);
        }

        // set actions
        for( size_t i = 0; i < numeric_features.size(); ++i ) {
            for( int t = 0; t <= d; ++t ) {
                std::string name = std::string("Set(") + numeric_features[i]->name() + "," + std::to_string(t) + ")";
                Action *action = new Action(name);
                for( size_t j = 0; j < numeric_features.size(); ++j )
                    action->add_precondition(q_features[j], false);
                action->add_precondition(bit_features[t], true);
                for( int j = t - 1; j >= 0; --j )
                    action->add_precondition(bit_features[j], false);
                action->add_effect(q_features[i], true);
                action->add_effect(bit_features[t], false);
                for( int j = t - 1; j >= 0; --j )
                    action->add_effect(bit_features[j], true);
                fond->add_action(action);
            }
        }

        // unset actions
        for( size_t i = 0; i < numeric_features.size(); ++i ) {
            std::string name = std::string("Unset(") + numeric_features[i]->name() + ")";
            Action *action = new Action(name);
            action->add_precondition(q_features[i], true);
            action->add_effect(q_features[i], false);
            fond->add_action(action);
        }

        // QNP actions
        for( size_t i = 0; i < actions_.size(); ++i ) {
            std::vector<const Feature*> features_increased;
            std::vector<const Feature*> features_decreased;
            for( size_t j = 0; j < actions_[i]->num_effects(); ++j ) {
                const Feature *feature = actions_[i]->effect(j).first;
                if( feature->numeric() && actions_[i]->effect(j).second )
                    features_increased.push_back(feature);
                else if( feature->numeric() )
                    features_decreased.push_back(feature);
            }

            if( !features_decreased.empty() ) {
                for( size_t j = 0; j < features_decreased.size(); ++j ) {
                    std::string name = actions_[i]->name();
                    if( features_decreased.size() > 1 )
                        name += std::string("_") + std::to_string(j);
                    Action *action = new Action(name);

                    // preconditions
                    for( size_t k = 0; k < actions_[i]->num_preconditions(); ++k ) {
                        const Feature *feature = actions_[i]->precondition(k).first;
                        bool value = actions_[i]->precondition(k).second;
                        action->add_precondition(fond->feature(feature->name()), value);
                    }
                    action->add_precondition(q_features[j], true);
                    for( size_t k = 0; k < features_increased.size(); ++k )
                        action->add_precondition(fond->feature(features_increased[k]->name()), false);

                    // effects
                    for( size_t k = 0; k < actions_[i]->num_effects(); ++k ) {
                        const Feature *feature = actions_[i]->effect(k).first;
                        bool value = actions_[i]->effect(k).second;
                        action->add_effect(fond->feature(feature->name()), value);
                    }
                    fond->add_action(action);
                }
            } else {
                Action *action = new Action(actions_[i]->name());

                // preconditions
                for( size_t k = 0; k < actions_[i]->num_preconditions(); ++k ) {
                    const Feature *feature = actions_[i]->precondition(k).first;
                    bool value = actions_[i]->precondition(k).second;
                    action->add_precondition(fond->feature(feature->name()), value);
                }

                // effects
                for( size_t k = 0; k < actions_[i]->num_effects(); ++k ) {
                    const Feature *feature = actions_[i]->effect(k).first;
                    bool value = actions_[i]->effect(k).second;
                    action->add_effect(fond->feature(feature->name()), value);
                }
                fond->add_action(action);
            }
        }

        return fond;
    }
#endif

    static Problem* read(std::istream &is) {
        std::string name;
        is >> name;
        Problem *qnp = new Problem(name);
        int num_features;
        is >> num_features;
        for( int i = 0; i < num_features; ++i ) {
            Feature *feature = Feature::read(is);
            qnp->add_feature(feature);
        }
        int num_actions;
        is >> num_actions;
        for( int i = 0; i < num_actions; ++i ) {
            Action *action = Action::read(is, qnp->feature_map());
            qnp->add_action(action);
        }
        return qnp;
    }

    void dump(std::ostream &os) const {
        os << name_ << std::endl << features_.size();
        for( size_t i = 0; i < features_.size(); ++i ) {
            os << " " << features_[i]->name() << " " << features_[i]->numeric();
        }
        os << std::endl << actions_.size() << std::endl;
        for( size_t i = 0; i < actions_.size(); ++i ) {
            os << *actions_[i];
        }
    }
    void PDDL_dump(std::ostream &os) const {
        os << "(define (domain " << PDDL_name(name_) << ")" << std::endl
           << "    (:requirements :non-deterministic)" << std::endl
           << "    (:types counter)" << std::endl
           << "    (:constants";

        for( size_t i = 0; i < features_.size(); ++i ) {
            if( features_[i]->numeric() )
                os << " " << PDDL_name(features_[i]->name());
        }
        os << " - counter)" << std::endl;

        os << "    (:predicates" << std::endl
           << "        (zero ?c - counter)" << std::endl
           << "        (q ?c - counter)" << std::endl;

        for( size_t i = 0; i < features_.size(); ++i ) {
            if( !features_[i]->numeric() && (features_[i]->PDDL_name().substr(0, 3) != "(q ") )
                os << "        " << features_[i]->PDDL_name() << std::endl;
        }
        os << "    )" << std::endl;

        if( !actions_.empty() ) os << std::endl;
        for( size_t i = 0; i < actions_.size(); ++i )
            actions_[i]->PDDL_dump(os);
        os << ")" << std::endl << std::endl;
    }
};

inline void Problem::clone(Problem *fond, const Action *action) const {
    Action *a = new Action(action->name());

    // preconditions
    for( size_t k = 0; k < action->num_preconditions(); ++k ) {
        const Feature *feature = action->precondition(k).first;
        bool value = action->precondition(k).second;
        a->add_precondition(fond->feature(feature->name()), value);
    }

    // effects
    for( size_t k = 0; k < action->num_effects(); ++k ) {
        const Feature *feature = action->effect(k).first;
        bool value = action->effect(k).second;
        a->add_effect(fond->feature(feature->name()), value);
    }
    fond->add_action(a);
}

inline void Problem::create_set_actions(Problem *fond,
                                        int d,
                                        int loop_nesting,
                                        std::vector<const Feature*> &numeric_features,
                                        std::vector<const Feature*> &boolean_features,
                                        std::vector<const Feature*> &bit_features,
                                        std::vector<const Feature*> &q_features) const {
    for( size_t i = 0; i < numeric_features.size(); ++i ) {
        for( int t = 0; t <= d; ++t ) {
            std::string name = std::string("Set(") + numeric_features[i]->name() + "," + std::to_string(t) + ")";
            Action *action = new Action(name);
            for( size_t j = 0; j < numeric_features.size(); ++j )
                action->add_precondition(q_features[j], false);
            action->add_precondition(bit_features[t], true);
            for( int j = t - 1; j >= 0; --j )
                action->add_precondition(bit_features[j], false);
            action->add_effect(q_features[i], true);
            action->add_effect(bit_features[t], false);
            for( int j = t - 1; j >= 0; --j )
                action->add_effect(bit_features[j], true);
            fond->add_action(action);
        }
    }
}

inline void Problem::translate(Problem *fond,
                               int d,
                               int loop_nesting,
                               std::vector<const Feature*> &numeric_features,
                               std::vector<const Feature*> &boolean_features,
                               std::vector<const Feature*> &bit_features,
                               std::vector<const Feature*> &q_features) const {
    if( loop_nesting == 0 ) {
        create_set_actions(fond, d, loop_nesting, numeric_features, boolean_features, bit_features, q_features);

        // unset actions
        for( size_t i = 0; i < numeric_features.size(); ++i ) {
            std::string name = std::string("Unset(") + numeric_features[i]->name() + ")";
            Action *action = new Action(name);
            action->add_precondition(q_features[i], true);
            action->add_effect(q_features[i], false);
            fond->add_action(action);
        }

        // QNP actions
        for( size_t i = 0; i < actions_.size(); ++i ) {
            const Action *action = actions_[i];

            std::vector<const Feature*> features_increased;
            std::vector<const Feature*> features_decreased;
            for( size_t j = 0; j < action->num_effects(); ++j ) {
                const Feature *feature = action->effect(j).first;
                if( feature->numeric() && action->effect(j).second )
                    features_increased.push_back(feature);
                else if( feature->numeric() )
                    features_decreased.push_back(feature);
            }

            if( !features_decreased.empty() ) {
                for( size_t j = 0; j < features_decreased.size(); ++j ) {
                    std::string name = action->name();
                    if( features_decreased.size() > 1 )
                        name += std::string("_") + std::to_string(j);
                    Action *a = new Action(name);

                    // preconditions
                    for( size_t k = 0; k < action->num_preconditions(); ++k ) {
                        const Feature *feature = action->precondition(k).first;
                        bool value = action->precondition(k).second;
                        a->add_precondition(fond->feature(feature->name()), value);
                    }
                    a->add_precondition(q_features[j], true);
                    for( size_t k = 0; k < features_increased.size(); ++k )
                        a->add_precondition(fond->feature(features_increased[k]->name()), false);

                    // effects
                    for( size_t k = 0; k < action->num_effects(); ++k ) {
                        const Feature *feature = action->effect(k).first;
                        bool value = action->effect(k).second;
                        a->add_effect(fond->feature(feature->name()), value);
                    }
                    fond->add_action(a);
                }
            } else {
                clone(fond, action);
            }
        }
    } else {
        std::cout << "**** Use loop_nesting=0" << std::endl;
        assert(0);
    }
}

inline Problem* Problem::translate(int d, int loop_nesting) const {
    // separate features
    std::vector<const Feature*> numeric_features;
    std::vector<const Feature*> boolean_features;
    for( size_t i = 0; i < features_.size(); ++i ) {
        if( features_[i]->numeric() )
            numeric_features.push_back(features_[i]);
        else
            boolean_features.push_back(features_[i]);
    }

    // create FOND
    Problem *fond = new Problem(std::string("FOND_") + name_ + "_" + std::to_string(d) + "_" + std::to_string(loop_nesting));
    for( size_t i = 0; i < numeric_features.size(); ++i )
        fond->add_feature(new Feature(numeric_features[i]->name(), true));
    for( size_t i = 0; i < boolean_features.size(); ++i )
        fond->add_feature(new Feature(boolean_features[i]->name(), false));

    // create new features
    std::vector<const Feature*> bit_features;
    for( int i = 0; i <= d; ++i ) {
        const Feature *feature = new Feature(std::string("bit(") + std::to_string(i) + ")", false);
        bit_features.push_back(feature);
        fond->add_feature(feature);
    }

    std::vector<const Feature*> q_features;
    for( size_t i = 0; i < numeric_features.size(); ++i ) {
        const Feature *feature = new Feature(PDDL_q(numeric_features[i]->name()), false, true);
        q_features.push_back(feature);
        fond->add_feature(feature);
    }

    // finish translation
    translate(fond, d, loop_nesting, numeric_features, boolean_features, bit_features, q_features);
    return fond;
}

}; // namespace QNP

inline std::ostream& operator<<(std::ostream &os, const QNP::Problem &qnp) {
    qnp.dump(os);
    return os;
}

#endif

