#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>

#include <blai/sample.h>
#include <blai/utils.h>
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

    // read arguments
    string prefix = argv[0];
    cout << "arguments: prefix=" << prefix << endl;

    // read data
    string matrix_filename = prefix + "/sat_matrix.dat";
    string transitions_filename = prefix + "/sat_transitions.dat";
    Sample::Sample sample(matrix_filename, transitions_filename);
    float preprocess_time = Utils::read_time_in_seconds() - start_time;

    // dump MaxSAT problem
    string wsat_filename = prefix + "_theory.wsat";
    CNFGenerator generator(sample);


    std::cout << Utils::blue() << "writing" << Utils::normal() << " '" << wsat_filename << "' ... " << std::endl;
    auto wsatstream = open_file(wsat_filename);
    // The generation of the CNF might already detect unsatisfiability
    bool unsat = generator.dump_weighted_max_sat_problem(wsatstream, true);
    wsatstream.close();

    if(unsat) {
        std::cout << Utils::warning() << "CNF theory is UNSAT" << std::endl;
    }

    float total_time = Utils::read_time_in_seconds() - start_time;

    cout << "stats:"
         << " preprocess-time " << preprocess_time
         << " total-time " << total_time
         << endl;

    return unsat ? 1 : 0;
}
