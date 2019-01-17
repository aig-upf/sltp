import logging
from collections import defaultdict
from tarski.dl import FeatureValueChange


class AbstractionValidator:
    """ Validate that certain abstractions are sound and complete wrt a given sample of state transitions """

    def __init__(self, model_cache, sample, state_ids):
        """ Build a validator from a given DLModelCache and sample of states. The states will be validated in the order
        given by `state_ids`, which is assumed to contain IDs of fully expanded states only """
        self.model_cache = model_cache
        self.sample = sample
        self.state_ids = state_ids

    def action_captures(self, models, s, sprime, action_effs, feature_set):
        """ Return true iff the abstract action captures (i.e. represents) transition (s, s') """
        for idx, f in feature_set:
            effect_of_action_on_f = action_effs[idx]
            x0 = self.get_possibly_cached_model(models, s).denotation(f)
            x1 = self.get_possibly_cached_model(models, sprime).denotation(f)
            diff = f.diff(x0, x1)
            if diff != effect_of_action_on_f:
                return False

        return True

    def action_set_captures(self, models, s, sprime, action_eff_set, feature_set):
        """ Return true iff some action in the given abstract action set captures transition (s, s') """
        return any(self.action_captures(models, s, sprime, effs, feature_set) for a, effs in action_eff_set.items())

    def get_possibly_cached_model(self, models, state_id):
        if state_id not in models:
            models[state_id] = self.model_cache.get_feature_model(state_id)
        return models[state_id]

    def find_flaws(self, abstraction, max_flaws):
        """ """
        sample = self.sample
        unsound, not_represented = set(), set()
        abstract_actions = abstraction["abstract_actions"]
        pool = abstraction["features"]

        feature_idx = self.compute_feature_idx(abstract_actions)
        models = dict()  # We will store here the model corresponding to each state

        # Map the idxs of selected features to the actual features
        selected_feature_idxs = abstraction["selected_features"]
        selected_features = [(f["idx"], abstraction["features"][f["idx"]]) for f in selected_feature_idxs]
        # assert all(str(f) == f2["name"] for f, f2 in zip(selected_features, selected_feature_idxs))

        logging.info("Abstract action set:\n{}".format(abstract_actions))

        assert all(s in sample.expanded for s in self.state_ids)  # state_ids assumed to contain only expanded states.
        for sid in self.state_ids:
            model = self.get_possibly_cached_model(models, sid)
            # Check soundness
            for action in abstract_actions:
                # We need to cast the actual feature value into a bool to compare it with the abstraction bool value
                is_applicable = all(bool(model.denotation(pool[idx])) is val for idx, val in action.preconditions)
                if is_applicable and \
                        not any(self.action_captures(models, sid, sprime, feature_idx[action], selected_features)
                                for sprime in sample.transitions[sid]):
                    # The abstract action is not sound
                    print("Action {} not *sound* wrt state {}".format(action, sid))
                    unsound.add(sid)

            # Check completeness
            if not all(self.action_set_captures(models, sid, sprime, feature_idx, selected_features)
                       for sprime in sample.transitions[sid]):
                # The abstract action is not complete
                print("Action set not *complete* wrt state {}".format(sid))
                not_represented.add(sid)

            if len(unsound) >= max_flaws:
                break

        return unsound if unsound else not_represented

    @staticmethod
    def compute_feature_idx(actions):
        index = dict()
        for a in actions:
            effects = [(eff.feature, eff.change) for eff in a.effects]
            index[a] = defaultdict(lambda: FeatureValueChange.NIL, effects)
        return index
