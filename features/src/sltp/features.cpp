
#include <sltp/features.hxx>

//#include <rapidcsv.h>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;


namespace SLTP { namespace DL {

        DL::State parse_state(unsigned id, const std::string &line) {
            std::vector<std::string> atoms;
            boost::split(atoms, line, boost::is_any_of("\t"));

            for (auto &atom_str:atoms) {
                std::vector<std::string> atom;
                boost::split(atom, atom_str, boost::is_any_of(","));  // Split e.g. "at,bob,shed" into a vector

//            for (const auto x:atom) std::cout << "Read atom: " << x << std::endl;
            }


            assert(0); // working on this
            return State(id);
        }

        std::vector<Object> parse_objects(const std::string &line) {
            std::vector<std::string> names;
            boost::split(names, line, boost::is_any_of(" \t"));

            unsigned object_id = 0;  // ATM we assign consecutive IDs to objects
            std::vector<Object> objects;
            for (auto &name:names) {
                objects.emplace_back(object_id, name);
                ++object_id;
            }

            return objects;
        }

        std::vector<Predicate> parse_predicates(const std::string &line) {
            std::vector<std::string> names;
            boost::split(names, line, boost::is_any_of(" \t"));

            unsigned id = 0;  // ATM we assign consecutive IDs
            std::vector<Predicate> predicates;
            for (auto &str:names) {
                std::vector<std::string> info;
                boost::split(info, str, boost::is_any_of("/"));  // Expected format: "clear/1"
                if (info.size() != 2) throw std::runtime_error("Wrong predicate format: " + str);

                predicates.emplace_back(id, info[0], (unsigned) std::stoi(info[1]));
                ++id;
            }

            return predicates;
        }


        Sample Sample::read(std::istream &is) {

            // First line: sample name or ID
            std::string sample_name;
            std::getline(is, sample_name); // Whole line is sample name

            // Second line: list of all predicate names
            std::string predicate_line;
            std::getline(is, predicate_line);
            std::vector<Predicate> predicates = parse_predicates(predicate_line);

            // Third line: list of all function names
            std::string function_line;
            std::getline(is, function_line);
            // TODO Not supporting functions yet - we don't do anything with this line

            // Fourth line: all possible object names in the sample (possibly from different problem instances)
            std::string object_line;
            std::getline(is, object_line);
            std::vector<Object> objects = parse_objects(object_line);


            std::vector<State> states;
            std::string line;
            for (unsigned id = 0; std::getline(is, line); ++id) {
                states.push_back(parse_state(id, line));
            }

            // TODO: Read ground predicates

            std::vector<const GroundedPredicate*> grounded_predicates;
            return Sample(sample_name, objects, predicates, grounded_predicates, states);
        }



    } };

