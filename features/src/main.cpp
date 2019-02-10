
#include <fstream>
#include <iostream>
#include <sltp/features.hxx>
#include <sltp/io.hxx>


using namespace std;

struct Options {
    std::string execname_;
    std::string worskspace_;
    unsigned complexity_bound_;

    Options(int argc, const char **argv) {
        // Print call
        cout << "Call:";
        for( int i = 0; i < argc; ++i )
            cout << " " << argv[i];
        cout << endl;

        execname_ = *argv;
        if( argc != 4 ) {
            print_usage(cout);
            exit(1);
        }

        // read arguments
        complexity_bound_ = (unsigned) stoul(argv[1], nullptr, 10);
        worskspace_ = string(argv[2]);
        cout << "Using workspace directory " << worskspace_ << endl;
    }

    void print_usage(ostream& os) const {
        os << endl
           << "Usage: " << execname_ << " <K> <input-file> <output-file>" << endl
           << endl
           << "where" << endl
           << "    <K> is the maximum feature complexity of the generated features"
           << "    <workspace> is the path to a directory with all the necessary input files "
              "                  from which to start the feature generation process." << endl
           << "    <output-file> is the path to the output CSV file where the denotation matrix of all generated features "
              "                  will be left." << endl
           << endl
           ;
    }
};

SLTP::DL::Sample parse_input_sample(const std::string& filename) {
    std::ifstream sample_file(filename);
    if (sample_file.fail()) throw std::runtime_error("Could not open filename '" + filename + "'");
    return SLTP::DL::Sample::read(sample_file);
}

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

    SLTP::DL::Sample sample = parse_input_sample(options.worskspace_ + "/sample.io");
    SLTP::DL::Factory factory("test", options.complexity_bound_);

    SLTP::DL::Cache cache;
    factory.generate_concepts(cache, &sample, true);
    factory.generate_features(cache, sample);

    return 0;
}

