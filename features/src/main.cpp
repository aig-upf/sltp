

#include <iostream>
#include <sltp/features.hxx>
//#include <sltp/io.hxx>

using namespace std;

struct Options {
    std::string execname_;
    std::string worskspace_;

    Options(int argc, const char **argv) {
        // Print call
        cout << "Call:";
        for( int i = 0; i < argc; ++i )
            cout << " " << argv[i];
        cout << endl;

        execname_ = *argv;
        if( argc != 3 ) {
            print_usage(cout);
            exit(1);
        }

        // read arguments
        worskspace_ = string(argv[1]);
        cout << "Using workspace directory " << worskspace_ << endl;
    }

    void print_usage(ostream& os) const {
        os << endl
           << "Usage: " << execname_ << " <input-file> <output-file>" << endl
           << endl
           << "where" << endl
           << "    <input-file> is the path to a CSV file containing the denotations of primitive concepts and roles "
              "                  from which to start the feature generation process." << endl
           << "    <output-file> is the path to the output CSV file where the denotation matrix of all generated features "
              "                  will be left." << endl
           << endl
           ;
    }
};

int main(int argc, const char **argv) {
    Options options(argc, argv);

    // Read primitive concepts and files from input file
    //auto primitive_denotations = sltp::io::read_denotation_matrix(options.worskspace + "/primitives.csv");
    //auto static_denotations = sltp::io::read_denotation_matrix(options.worskspace + "/statics.csv");

    // Generate concepts based on grammar, pruning duplicate concepts

    // Generate Features and prune duplicates

    // Generate output "state x feature" matrix
    // We will likely need to output somehow the structure of the (non-redundant) concepts and features,
    // so that they can be reconstructed from python.

    return 0;
}

