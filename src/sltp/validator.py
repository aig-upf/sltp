import logging
from collections import defaultdict

from tarski.dl import FeatureValueChange

from .util.tools import Abstraction


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
        for f in feature_set:
            effect_of_action_on_f = action_effs[f.id]
            x0 = f.denotation(self.get_possibly_cached_model(models, s))
            x1 = f.denotation(self.get_possibly_cached_model(models, sprime))
            diff = f.feature.diff(x0, x1)
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

    def is_applicable(self, models, sid, action):
        """ Return true if the given abstract action is applicable in the given state """
        model = self.get_possibly_cached_model(models, sid)
        # We need to cast the actual feature value into a bool to compare it with the abstraction bool value
        return all(bool(feat.denotation(model)) is val for feat, val in action.preconditions)

    def find_flaws(self, abstraction, max_flaws, check_completeness=True):
        """ """
        assert isinstance(abstraction, Abstraction)
        sample = self.sample
        unsound, not_represented = set(), set()
        feature_idx = self.compute_feature_idx(abstraction.actions)
        models = dict()  # We will store here the model corresponding to each state

        assert all(s in sample.expanded for s in self.state_ids)  # state_ids assumed to contain only expanded states.
        for sid in self.state_ids:

            # Check soundness
            for action in abstraction.actions:
                if self.is_applicable(models, sid, action) and \
                        not any(self.action_captures(models, sid, sprime, feature_idx[action], abstraction.features)
                                for sprime in sample.transitions[sid]):
                    # The abstract action is not sound
                    logging.info("Action {} *unsound* wrt state #{}".format(action, sid))
                    logging.info("State #{} is: {}".format(sid, sample.states[sid]))
                    children = map(str, (sample.states[sprime] for sprime in sample.transitions[sid]))
                    logging.info("Children of #{} are:\n\t{}".format(sid, "\n\t".join(children)))

                    unsound.add(sid)
                    break

            # Check completeness
            if check_completeness:
                if not all(self.action_set_captures(models, sid, sprime, feature_idx, abstraction.features)
                           for sprime in sample.transitions[sid]):
                    # The abstract action is not complete
                    logging.debug("Action set not *complete* wrt state {}".format(sid))
                    not_represented.add(sid)

            if len(unsound) >= max_flaws:
                break

        return unsound if unsound else not_represented

    @staticmethod
    def compute_feature_idx(actions):
        """ Return a mapping from abstract actions to the qualitative effect (ADD, DEL, INC, DEC, NONE)
        that each effect of the action has on each feature index.
        """
        index = dict()
        for a in actions:
            effects = [(eff.feature.id, eff.change) for eff in a.effects]
            index[a] = defaultdict(lambda: FeatureValueChange.NIL, effects)
        return index
