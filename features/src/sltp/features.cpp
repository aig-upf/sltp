
#include <sltp/features.hxx>

//#include <rapidcsv.h>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;


namespace SLTP { namespace DL {

        DL::State parse_state(unsigned id, const std::string &state_line) {
            std::vector<std::string> atoms;
            boost::split(atoms, state_line, boost::is_any_of("\t"));

            for (auto &atom_str:atoms) {
                std::vector<std::string> atom;
                boost::split(atom, atom_str, boost::is_any_of(","));  // Split e.g. "at,bob,shed" into a vector

//            for (const auto x:atom) std::cout << "Read atom: " << x << std::endl;
            }



        }

        std::vector<Object> parse_objects(const std::string &object_line) {
            std::vector<std::string> names;
            boost::split(names, object_line, boost::is_any_of(" \t"));

            std::vector<Object> objects;
            for (auto &name:names) {
                objects.emplace_back(name);
            }

            return objects;
        }


        Sample Sample::read(std::istream &is) {



            // First line: sample name or ID
            std::string sample_name;
            std::getline(is, sample_name);

            // Second line: all possible object names in the sample (possibly from different problem instances)
            std::string object_line;
            std::getline(is, object_line);

            std::vector<Object> objects = parse_objects(object_line);



//            for (; ; ) {
//                DL::State s = parse_state(id, line);
//
//
//                std::istringstream iss(line);
//                int a, b;
//                if (!(iss >> a >> b)) { break; } // error
//
//                // process pair (a,b)
//
//
//            }

            std::vector<State> states;
            std::string line;
            for (unsigned id = 0; std::getline(is, line); ++id) {
                states.push_back(parse_state(id, line));
            }


//            throw std::runtime_error("Work in Progress");

            std::vector<const Predicate*> predicates;
            std::vector<const GroundedPredicate*> grounded_predicates;
            return Sample(sample_name, objects, predicates, grounded_predicates, states);
        }



    } };

