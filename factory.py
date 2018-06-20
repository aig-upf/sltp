
import tarski as tsk
import logging

from tarski.syntax import builtins

from tarski.dl import Concept, Role, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, RestrictRole, AndConcept, EmptyConcept, CompositionRole, SingletonConcept


def filter_subnodes(elem, t):
    return list(filter(lambda x: type(x) == t, elem.flatten()))


class SyntacticFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage):
        self.language = language
        self.universe_sort = language.get_sort('object')
        self.top = UniversalConcept(self.universe_sort.name)
        self.bot = EmptyConcept(self.universe_sort.name)

    def generate_primitives_from_language(self):
        assert len(self.language.functions) == 0

        concepts, roles, primitive_atoms = [], [], []
        for predicate in self.language.predicates:
            if builtins.is_builtin_predicate(predicate):
                # Skip "=" and other built-in symbols
                # TODO We might want to revise this and allow for certain builtins
                continue

            name = predicate.symbol

            if predicate.arity == 0:
                primitive_atoms.append(name)
            elif predicate.arity == 1:
                concepts.append(BasicConcept(predicate))
            elif predicate.arity == 2:
                roles.append(BasicRole(predicate))
            else:
                print("Predicate {} with arity > 2 ignored".format(predicate))

        for c in self.language.constants():
            concepts.append(SingletonConcept(c.symbol, c.sort))

        logging.info('Primitive (nullary) atoms : {}'.format(", ".join(map(str, primitive_atoms))))
        logging.info('Primitive (unary) concepts: {}'.format(", ".join(map(str, concepts))))
        logging.info('Primitive (binary) roles  : {}'.format(", ".join(map(str, roles))))

        # TODO At the moment we do nothing with nullary atoms (!!!)
        # if primitive_atoms:
        #     logging.critical("Support for primitive nullary atoms not yet implemented!!")

        return concepts, roles, primitive_atoms

    def create_exists_concept(self, role: Role, concept: Concept):
        # exists(R.C) = { x | exists y R(x,y) and C(y) }

        result = ExistsConcept(role, concept)
        s1, s2 = role.sort

        if concept == self.bot:
            logging.debug('Concept "{}" is statically empty'.format(result))
            return None

        # TODO ADD: If C is a sort-concept of the same sort than s2, then the concept will be equiv to exist(R.True)
        if not self.language.are_vertically_related(s2, concept.sort):
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        if isinstance(role, RestrictRole) and concept == self.top:
            return None  # Is syntactically equivalent to a simpler exists concept

        return result

    def create_not_concept(self, concept: Concept):
        if isinstance(concept, NotConcept):
            return None

        result = NotConcept(concept, self.universe_sort)
        return result

    def create_forall_concept(self, role: Role, concept: Concept):
        # forall(R.C) = { x | forall y R(x,y) implies C(y) }

        result = ForallConcept(role, concept)
        s1, s2 = role.sort

        if isinstance(concept, UniversalConcept):
            # logging.debug('Concept "{}" equivalent to simpler concept'.format(result))
            return None

        if not self.language.are_vertically_related(s2, concept.sort):
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_and_concept(self, c1: Concept, c2: Concept):
        # C1 AND C2 = { x | C1(x) AND C2(x) }

        sort = self.language.most_restricted_type(c1.sort, c2.sort)

        if c1 == c2:
            return None  # No sense in C and C

        if c1 in (self.top, self.bot) or c2 in (self.top, self.bot):
            logging.debug('AND of {} and {} pruned, no sense in AND\'ing with top or bot'.format(c1, c2))
            return None

        if sort is None:
            # i.e. c1 and c2 are disjoint types
            logging.debug('AND of {} and {} pruned for type-inconsistency reasons'.format(c1, c2))
            return None

        return AndConcept(c1, c2, sort)

    def create_equal_concept(self, r1: Role, r2: Role):
        assert isinstance(r1, Role) and isinstance(r2, Role)
        # The sort of the resulting concept will be the most restricted sort between the left sorts of the two roles
        sort = self.language.most_restricted_type(r1.sort[0], r2.sort[0])
        result = EqualConcept(r1, r2, sort)

        if not self.language.are_vertically_related(r1.sort[0], r2.sort[0]) or \
                not self.language.are_vertically_related(r1.sort[1], r2.sort[1]):
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_restrict_role(self, r: Role, c: Concept):

        result = RestrictRole(r, c)
        if not self.language.are_vertically_related(r.sort[1], c.sort):
            logging.debug('Role "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        if isinstance(c, UniversalConcept) or c == self.bot:
            logging.debug('Role "{}" pruned; no sense in restricting to top / bot concepts'.format(result))
            return None

        if isinstance(r, RestrictRole):
            logging.debug('Role "{}" pruned; no direct nesting of restrictions'.format(result))
            return None

        return result

    def create_composition_role(self, r1: Role, r2: Role):

        # Compose only on primitives or their inversions
        # if (not isinstance(r1, (BasicRole, InverseRole)) or
        #    not isinstance(r2, (BasicRole, InverseRole))):
        #     return None

        result = CompositionRole(r1, r2)

        if not self.language.are_vertically_related(r1.sort[1], r2.sort[0]):
            logging.debug('Role "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        num_comp = len(filter_subnodes(result, CompositionRole))
        if num_comp > 2:
            logging.debug('Role "{}" pruned: number of compositions ({}) exceeds threshold'.format(result, num_comp))
            return None

        return result
