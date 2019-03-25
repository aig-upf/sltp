#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>

#include <blai/sample.h>
#include <blai/utils.h>
#include <blai/matrix.h>
#include <blai/transitions.h>
#include <cnf/generator.h>
#include <common/utils.h>

using namespace std;

void usage(ostream &os, const string &name) {
    cout << endl
         << "usage: " << name << " <prefix>" << endl
         << endl
         << "where" << endl
         << "    <prefix> is prefix for all files" << endl
         << endl
         << "Constructs a Weighted MaxSAT problem from matrix and transition files." << endl
         << "The output is a .wsat file that can be fed directly into the solver." << endl
         << endl;
}

int main(int argc, const char **argv) {
    // print call and start clock
    cout << "call: " << Utils::cmdline(argc, argv) << endl;
    float start_time = Utils::read_time_in_seconds();

    // read executable name
    string name(*argv);

    // read options
    for( ++argv, --argc; (argc > 0) && (**argv == '-'); ++argv, --argc ) {
        if( string(*argv) == "--" ) {
            ++argv;
            --argc;
            break;
        } else {
            cout << Utils::error() << "unrecognized option '" << *argv << "'" << endl;
            usage(cout, name);
            exit(0);
        }
    }

    // check we have enough arguments
    if( argc < 1 ) {
        cout << Utils::error() << "insufficient arguments" << endl;
        usage(cout, name);
        exit(0);
    }

    bool verbose = true;

    // read arguments
    string prefix = argv[0];
    cout << "arguments: prefix=" << prefix << endl;

    // Read feature matrix
    std::string matrix_filename = prefix + "/feature-matrix.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << matrix_filename << std::endl;
    auto ifs_matrix = get_ifstream(matrix_filename);
    auto matrix = Sample::Matrix::read_dump(ifs_matrix, verbose);
    ifs_matrix.close();

    // Read transition data
    std::string transitions_filename = prefix + "/sat_transitions.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << transitions_filename << std::endl;
    auto ifs_transitions = get_ifstream(transitions_filename);
    auto transitions = Sample::Transitions::read_dump(ifs_transitions, verbose);
    ifs_transitions.close();


    Sample::Sample sample(std::move(matrix), std::move(transitions), verbose);
    float preprocess_time = Utils::read_time_in_seconds() - start_time;

    // dump MaxSAT problem
    string top_filename = prefix + "/top.dat";
    string wsat_filename = prefix + "/theory.wsat";
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
