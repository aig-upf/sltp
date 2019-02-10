
#include <sltp/features.hxx>

//#include <rapidcsv.h>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;


namespace SLTP { namespace DL {

        DL::GroundedPredicate
        parse_atom(const std::string &line, const Sample::ObjectIndex& oidx, const Sample::PredicateIndex& pidx) {
            std::vector<std::string> atom;
            boost::split(atom, line, boost::is_any_of(","));  // Split e.g. "at,bob,shed" into a vector
            if (atom.empty()) throw std::runtime_error("Wrong atom format: " + line);


            std::vector<object_id> objects;
            objects.reserve(atom.size()-1);
            for (unsigned i = 1; i < atom.size(); ++i) { // Start iterating at 1 to skip the predicate name
                objects.push_back(oidx.at(atom[i]));
            }

            return GroundedPredicate(pidx.at(atom[0]), objects);
        }

        std::vector<DL::GroundedPredicate>
        parse_atoms(const std::string &line,
                const Sample::ObjectIndex& oidx, const Sample::PredicateIndex& pidx, Sample::AtomIndex& aidx) {
            std::vector<std::string> atomlist;
            boost::split(atomlist, line, boost::is_any_of("\t"));

            std::vector<DL::GroundedPredicate> atoms;
            atoms.reserve(atomlist.size());

            for (auto &atom_str:atomlist) {
                auto aid = (atom_id) atoms.size();
                atoms.push_back(parse_atom(atom_str, oidx, pidx));

                // Update the atom index
                aidx.insert(std::make_pair(atoms.back().data(), aid));
            }

            return atoms;
        }

        DL::State parse_state(unsigned id, const std::string &line,
                const Sample::ObjectIndex& oidx, const Sample::PredicateIndex& pidx, Sample::AtomIndex& aidx) {

            std::vector<std::string> atomlist;
            boost::split(atomlist, line, boost::is_any_of("\t"));

            std::vector<atom_id> atom_ids;
            atom_ids.reserve(atomlist.size());


            for (auto &atom_str:atomlist) {
                auto atom = parse_atom(atom_str, oidx, pidx);
                atom_ids.push_back(aidx.at(atom.data()));
            }

            return State(id, atom_ids);
        }

        std::vector<Object> parse_objects(const std::string &line, Sample::ObjectIndex& oidx) {
            std::vector<std::string> names;
            boost::split(names, line, boost::is_any_of(" \t"));

            std::vector<Object> objects;
            for (auto &name:names) {
                auto oid = (object_id) oidx.size();
                objects.emplace_back(oid, name);
                oidx.insert(std::make_pair(name, oid));
            }

            return objects;
        }

        std::vector<Predicate> parse_predicates(const std::string &line, Sample::PredicateIndex& pidx) {
            std::vector<std::string> names;
            boost::split(names, line, boost::is_any_of(" \t"));

            std::vector<Predicate> predicates;
            for (auto &str:names) {
                std::vector<std::string> info;
                boost::split(info, str, boost::is_any_of("/"));  // Expected format: "clear/1"
                if (info.size() != 2) throw std::runtime_error("Wrong predicate format: " + str);

                auto pid= (object_id) pidx.size();
                predicates.emplace_back(pid, info[0], (unsigned) std::stoi(info[1]));
                pidx.insert(std::make_pair(info[0], pid));

            }

            return predicates;
        }


        Sample Sample::read(std::istream &is) {
            ObjectIndex oidx;
            PredicateIndex pidx;
            AtomIndex aidx;

            // First line: sample name or ID
            std::string sample_name;
            std::getline(is, sample_name); // Whole line is sample name

            // Second line: list of all predicate names
            std::string predicate_line;
            std::getline(is, predicate_line);
            std::vector<Predicate> predicates = parse_predicates(predicate_line, pidx);

            // Third line: list of all function names
            std::string function_line;
            std::getline(is, function_line);
            // Note: Not supporting functions yet - we don't do anything with this line

            // Fourth line: all possible object names in the sample (possibly from different problem instances)
            std::string object_line;
            std::getline(is, object_line);
            std::vector<Object> objects = parse_objects(object_line, oidx);

            // Fifth line: all possible atoms, i.e. grounded predicates (possibly from different problem instances)
            std::string atom_line;
            std::getline(is, atom_line);
            std::vector<GroundedPredicate> atoms = parse_atoms(atom_line, oidx, pidx, aidx);

            // Remaining lines: one line per state with all of the atoms in the state
            std::vector<State> states;
            std::string line;
            for (unsigned id = 0; std::getline(is, line); ++id) {
                states.push_back(parse_state(id, line, oidx, pidx, aidx));
            }

            return Sample(sample_name, objects, predicates, atoms, states, oidx, pidx, aidx);
        }


} }; // namespaces

