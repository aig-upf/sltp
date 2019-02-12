
#include <sltp/features.hxx>

//#include <rapidcsv.h>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;


namespace SLTP { namespace DL {

    DL::GroundedPredicate parse_atom(const std::string &line,
                                     const Sample::ObjectIndex &object_index,
                                     const Sample::PredicateIndex &predicate_index) {
        std::vector<std::string> atom;
        boost::split(atom, line, boost::is_any_of(","));  // Split e.g. "at,bob,shed" into a vector
        if( atom.empty() ) throw std::runtime_error("Wrong atom format: " + line);

        std::vector<object_id> objects;
        objects.reserve(atom.size() - 1);
        for( int i = 1; i < atom.size(); ++i ) { // Start iterating at 1 to skip the predicate name
            objects.push_back(object_index.at(atom[i]));
        }
        return GroundedPredicate(predicate_index.at(atom[0]), std::move(objects));
    }

    std::vector<DL::GroundedPredicate> parse_atoms(const std::string &line,
                                                   const Sample::ObjectIndex &object_index,
                                                   const Sample::PredicateIndex &predicate_index,
                                                   Sample::AtomIndex &atom_index) {
        std::vector<std::string> atom_list;
        boost::split(atom_list, line, boost::is_any_of("\t"));

        std::vector<DL::GroundedPredicate> atoms;
        atoms.reserve(atom_list.size());

        for( auto &atom_str : atom_list ) {
            auto aid = (atom_id) atoms.size();
            atoms.emplace_back(std::move(parse_atom(atom_str, object_index, predicate_index)));

            // Update the atom index
            atom_index.insert(std::make_pair(atoms.back().data(), aid));
        }

        return atoms;
    }

    DL::State parse_state(const Instance &ins,
                          unsigned id,
                          const std::string &line,
                          const Sample::ObjectIndex &object_index,
                          const Sample::PredicateIndex &predicate_index,
                          Sample::AtomIndex &atom_index) {
        std::vector<std::string> atom_list;
        boost::split(atom_list, line, boost::is_any_of("\t"));

        std::vector<atom_id> atom_ids;
        atom_ids.reserve(atom_list.size());

        for( auto &atom_str : atom_list ) {
            auto atom = std::move(parse_atom(atom_str, object_index, predicate_index));
            atom_ids.push_back(atom_index.at(atom.data()));
        }

        return State(ins, id, std::move(atom_ids));
    }

    std::vector<Object> parse_objects(const std::string &line,
                                      Sample::ObjectIndex &object_index) {
        std::vector<std::string> names;
        boost::split(names, line, boost::is_any_of(" \t"));

        std::vector<Object> objects;
        for( auto &name : names ) {
            auto oid = (object_id) object_index.size();
            objects.emplace_back(oid, name);
            object_index.insert(std::make_pair(name, oid));
        }

        return objects;
    }

    std::vector<Predicate> parse_predicates(const std::string &line,
                                            Sample::PredicateIndex& predicate_index) {
        std::vector<std::string> names;
        boost::split(names, line, boost::is_any_of(" \t"));

        std::vector<Predicate> predicates;
        for( auto &str : names ) {
            std::vector<std::string> info;
            boost::split(info, str, boost::is_any_of("/"));  // Expected format: "clear/1"
            if( info.size() != 2 ) throw std::runtime_error("Wrong predicate format: " + str);

            auto pid = (object_id) predicate_index.size();
            predicates.emplace_back(pid, info[0], (unsigned) std::stoi(info[1]));
            predicate_index.insert(std::make_pair(info[0], pid));
        }

        return predicates;
    }

    const Sample* Sample::read(std::istream &is) {
        ObjectIndex object_index;
        PredicateIndex predicate_index;
        AtomIndex atom_index;

        // First line: sample name or ID
        std::string sample_name;
        std::getline(is, sample_name); // Whole line is sample name

        // Second line: list of all predicate names
        std::string predicate_line;
        std::getline(is, predicate_line);
        std::vector<Predicate> predicates(std::move(parse_predicates(predicate_line, predicate_index)));

        // Third line: list of all function names
        std::string function_line;
        std::getline(is, function_line);
        // Note: Not supporting functions yet - we don't do anything with this line

        // Fourth line: all possible object names in the sample (possibly from different problem instances)
        std::string object_line;
        std::getline(is, object_line);
        std::vector<Object> objects(std::move(parse_objects(object_line, object_index)));

        // Fifth line: all possible atoms, i.e. grounded predicates (possibly from different problem instances)
        std::string atom_line;
        std::getline(is, atom_line);
        std::vector<GroundedPredicate> atoms(std::move(parse_atoms(atom_line, object_index, predicate_index, atom_index)));

        // Remaining lines: one line per state with all of the atoms in the state
        Instance ins;
        std::vector<State> states;
        std::string line;
        for( unsigned id = 0; std::getline(is, line); ++id )
            states.emplace_back(std::move(parse_state(ins, id, line, object_index, predicate_index, atom_index)));

        return new Sample(sample_name,
                          std::move(objects),
                          std::move(predicates),
                          std::move(atoms),
                          std::move(states),
                          std::move(object_index),
                          std::move(predicate_index),
                          std::move(atom_index));
    }

} }; // namespaces

