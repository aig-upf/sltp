
#include <sltp/features.hxx>

//#include <rapidcsv.h>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;


namespace SLTP { namespace DL {

    DL::GroundedPredicate parse_atom(const std::string &line,
                                     const Instance::ObjectIndex &object_index,
                                     const Sample::PredicateIndex &predicate_index) {
        std::vector<std::string> atom;
        boost::split(atom, line, boost::is_any_of(","));  // Split e.g. "at,bob,shed" into a vector
        if( atom.empty() ) throw std::runtime_error("Wrong atom format: " + line);

        std::vector<object_id_t> objects;
        objects.reserve(atom.size() - 1);
        for( int i = 1; i < atom.size(); ++i ) { // Start iterating at 1 to skip the predicate name
            objects.push_back(object_index.at(atom[i]));
        }
        return GroundedPredicate(predicate_index.at(atom[0]), std::move(objects));
    }

    std::vector<DL::GroundedPredicate> parse_atoms(const std::string &line,
                                                   const Instance::ObjectIndex &object_index,
                                                   const Sample::PredicateIndex &predicate_index,
                                                   Instance::AtomIndex &atom_index) {
        std::vector<std::string> atom_list;
        boost::split(atom_list, line, boost::is_any_of("\t"));

        std::vector<DL::GroundedPredicate> atoms;
        atoms.reserve(atom_list.size());

        for( auto &atom_str : atom_list ) {
            auto aid = (atom_id_t) atoms.size();
            atoms.emplace_back(std::move(parse_atom(atom_str, object_index, predicate_index)));

            // Update the atom index
            atom_index.insert(std::make_pair(atoms.back().data(), aid));
        }

        return atoms;
    }

    DL::State parse_state(const Instance &ins,
                          unsigned id,
                          const std::string &line,
                          const Instance::ObjectIndex &object_index,
                          const Sample::PredicateIndex &predicate_index,
                          Instance::AtomIndex &atom_index) {
        std::vector<std::string> atom_list;
        boost::split(atom_list, line, boost::is_any_of("\t"));

        std::vector<atom_id_t> atom_ids;
        atom_ids.reserve(atom_list.size());

        for( auto &atom_str : atom_list ) {
            auto atom = std::move(parse_atom(atom_str, object_index, predicate_index));
            atom_ids.push_back(atom_index.at(atom.data()));
        }

        return State(ins, id, std::move(atom_ids));
    }

    std::vector<Object> parse_objects(const std::string &line, Instance::ObjectIndex &object_index) {
        std::vector<std::string> names;
        boost::split(names, line, boost::is_any_of(" \t"));

        std::vector<Object> objects;
        for( auto &name : names ) {
            auto oid = (object_id_t) object_index.size();
            objects.emplace_back(oid, name);
            object_index.insert(std::make_pair(name, oid));
        }

        return objects;
    }

    std::vector<Predicate> parse_predicates(const std::string &line, Sample::PredicateIndex& predicate_index) {
        std::vector<std::string> names;
        boost::split(names, line, boost::is_any_of(" \t"));

        std::vector<Predicate> predicates;
        for( auto &str : names ) {
            std::vector<std::string> info;
            boost::split(info, str, boost::is_any_of("/"));  // Expected format: "clear/1"
            if( info.size() != 2 ) throw std::runtime_error("Wrong predicate format: " + str);

            auto pid = (object_id_t) predicate_index.size();
            predicates.emplace_back(pid, info[0], (unsigned) std::stoi(info[1]));
            predicate_index.insert(std::make_pair(info[0], pid));
        }

        return predicates;
    }

    const Sample* Sample::read(std::istream &is) {
        //ObjectIndex object_index;
        PredicateIndex predicate_index;
        //AtomIndex atom_index;

        // First line: sample name or ID
        std::string sample_name;
        std::getline(is, sample_name); // Whole line is sample name

        // Second line: list of all predicate names (shared by all instances)
        std::string predicate_line;
        std::getline(is, predicate_line);
        std::vector<Predicate> predicates(std::move(parse_predicates(predicate_line, predicate_index)));

        // Third line: list of all function names (shared by all instances)
        std::string function_line;
        std::getline(is, function_line);
        // Note: Not supporting functions yet - we don't do anything with this line

        // Fourth line: number of instances in file
        std::string num_instances_line;
        std::getline(is, num_instances_line);
        int num_instances = atoi(num_instances_line.c_str());

        // the rest of the file contains instances. Each instance
        // consists of a header containing instance name and number
        // of states in the instance, a list of object names, a
        // list of grounded predicates, and the list of states (one
        // per line)
        std::vector<Instance> instances;
        std::vector<State> states;
        for( int i = 0; i < num_instances; ++i ) {
            std::cout << "hola" << std::endl;
            // First line in instance: header with name and number of states
            std::string instance_header_line;
            std::getline(is, instance_header_line);
            std::vector<std::string> tokens;
            boost::split(tokens, instance_header_line, boost::is_any_of(" "));
            if( tokens.size() != 2 ) throw std::runtime_error("Wrong instance header format: " + instance_header_line);

            std::string instance_name = tokens[0];
            int num_states = atoi(tokens[1].c_str());

            // indices
            Instance::ObjectIndex object_index;
            Instance::AtomIndex atom_index;

            // Second line in instance: all possible object names in this instance
            std::string object_line;
            std::getline(is, object_line);
            std::vector<Object> objects(std::move(parse_objects(object_line, object_index)));

            // Third line in instance: all possible atoms in instance, i.e. grounded predicates
            std::string atom_line;
            std::getline(is, atom_line);
            std::vector<GroundedPredicate> atoms(std::move(parse_atoms(atom_line, object_index, predicate_index, atom_index)));

            // create instance
            Instance instance(instance_name, std::move(objects), std::move(atoms), std::move(object_index), std::move(atom_index));

            // Remaining lines in instance: one line per state with all of the atoms in the state
            int base_id = states.size();
            for( int id = 0; id < num_states; ++id ) {
                std::string state_line;
                std::getline(is, state_line);
                states.emplace_back(std::move(parse_state(instance, base_id + id, state_line, object_index, predicate_index, atom_index)));
            }

            // insert instance into list of instances
            instances.emplace_back(std::move(instance));
        }

        // create and return sample
        const Sample *sample = new Sample(sample_name,
                                          std::move(predicates),
                                          std::move(instances),
                                          std::move(states),
                                          std::move(predicate_index));
        return sample;
    }

} }; // namespaces

