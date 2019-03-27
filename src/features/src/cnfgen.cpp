#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <blai/sample.h>
#include <blai/utils.h>
#include <blai/matrix.h>
#include <blai/transitions.h>
#include <cnf/generator.h>
#include <common/utils.h>


using namespace std;
namespace po = boost::program_options;

struct Options {
    std::string workspace;
    bool prune_redundant_states;
    bool verbose;

    Options(int argc, const char **argv) {
        po::options_description description("Generate a weighted max-sat instance from given feature "
                                            "matrix and transition samples.\n\nOptions");

        description.add_options()
                ("help,h", "Display this help message and exit.")
                ("verbose,v", "Show extra debugging messages.")

                ("workspace,w", po::value<std::string>()->required(),
                        "Directory where the input files (feature matrix, transition sample) reside, "
                        "and where the output .wsat file will be left.")

                ("prune-redundant-states,p",
                        "Whether to prune those states that appear redundant for the given feature pool.")
                ;

        po::variables_map vm;

        try {
            po::store(po::command_line_parser(argc, argv).options(description).run(), vm);

            if (vm.count("help")) {
                std::cout << description << "\n";
                exit(0);
            }
            po::notify(vm);
        } catch (const std::exception &ex) {
            std::cout << "Error with command-line options:" << ex.what() << std::endl;
            std::cout << std::endl << description << std::endl;
            exit(1);
        }

        workspace = vm["workspace"].as<std::string>();
        prune_redundant_states = vm.count("prune-redundant-states") > 0;
        verbose = vm.count("verbose") > 0;
    }
};


int main(int argc, const char **argv) {
    float start_time = Utils::read_time_in_seconds();
    Options options(argc, argv);

    // Read feature matrix
    std::string matrix_filename = options.workspace + "/feature-matrix.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << matrix_filename << std::endl;
    auto ifs_matrix = get_ifstream(matrix_filename);
    auto matrix = Sample::Matrix::read_dump(ifs_matrix, options.verbose);
    ifs_matrix.close();

    // Read transition data
    std::string transitions_filename = options.workspace + "/sat_transitions.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << transitions_filename << std::endl;
    auto ifs_transitions = get_ifstream(transitions_filename);
    auto transitions = Sample::Transitions::read_dump(ifs_transitions, options.verbose);
    ifs_transitions.close();


    Sample::Sample sample(std::move(matrix), std::move(transitions), options.verbose);
    float preprocess_time = Utils::read_time_in_seconds() - start_time;

    // dump MaxSAT problem
    string top_filename = options.workspace + "/top.dat";
    string wsat_filename = options.workspace + "/theory.wsat";
    string wsat_filename_tmp = wsat_filename + ".tmp";
    CNFGenerator gen(sample);

    auto wsatstream = get_ofstream(wsat_filename_tmp);
    auto res = gen.write_maxsat(wsatstream, true);
    wsatstream.close();

    // return {true, writer.top(), writer.nvars(), writer.nclauses()};
    const auto& writer = res.second;

    auto topstream = get_ofstream(top_filename);
    topstream << writer.top() << " "
              << writer.nvars() << " "
              << writer.nclauses();
    topstream.close();

    bool unsat = std::get<0>(res);
    if(unsat) {
        // The generation of the CNF might have already detected unsatisfiability
        std::cout << Utils::warning() << "CNF theory is UNSAT" << std::endl;
    }

    float total_time = Utils::read_time_in_seconds() - start_time;

    std::cout << "Runtime stats:"  << " preprocess-time " << preprocess_time << " total-time " << total_time << std::endl;
    std::cout << "CNF Theory: "  << writer.nvars() << " vars + " << writer.nclauses() << " clauses" << std::endl;
    std::cout << "Clause breakdown: " << std::endl;
    std::cout << "Selected(f):" <<  gen.n_selected_clauses << std::endl;
    std::cout << "D2-definition:" <<  gen.n_d2_clauses << std::endl;
    std::cout << "Bridge:" <<  gen.n_bridge_clauses << std::endl;
    std::cout << "Goal-distinguishing:" <<  gen.n_goal_clauses << std::endl;
    std::cout << "TOTAL:" <<  gen.n_selected_clauses + gen.n_d2_clauses + gen.n_bridge_clauses + gen.n_goal_clauses << std::endl;

    return unsat ? 1 : 0;
}
