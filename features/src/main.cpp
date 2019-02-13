
#include <iostream>
#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>

#include <sltp/features.hxx>

using namespace std;

struct Options {
    string execname_;
    string workspace_;
    unsigned complexity_bound_;

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
        complexity_bound_ = (unsigned) stoul(argv[1], nullptr, 10);
        workspace_ = string(argv[2]);
        //cout << "Using workspace directory '" << workspace_ << "' and output file '" << output_file_ << "'" << endl;
    }

    void print_usage(ostream &os) const {
        os << endl
           << "Usage: " << execname_ << " <K> <input-file> <output-file>" << endl
           << endl
           << "where" << endl
           << "    <K> is the maximum feature complexity of the generated features."
           << endl
           << "    <workspace> is the path to a directory with all the necessary input"
           << endl
           << "                files from which to start the feature generation process."
           << endl
           ;
    }
};

SLTP::DL::Sample parse_input_sample(const string &filename) {
    ifstream sample_file(filename);
    if( sample_file.fail() ) throw runtime_error("Could not open filename '" + filename + "'");
    return SLTP::DL::Sample::read(sample_file);
}

std::vector<std::string> parse_nominals(const string &filename) {
    ifstream file(filename);
    if( file.fail() ) throw runtime_error("Could not open filename '" + filename + "'");

    std::string nominals_line;
    std::getline(file, nominals_line);
    if (nominals_line.empty()) return {};
    std::vector<std::string> nominals;
    boost::split(nominals, nominals_line, boost::is_any_of(" \t"));
    return nominals;
}

void output_results(const Options &options,
                    const SLTP::DL::Factory &factory,
                    const SLTP::DL::Sample &sample,
                    const SLTP::DL::Cache &cache) {

    // Print feature matrix
    string output(options.workspace_ + "/feature-matrix.io");
    ofstream output_file(output);
    if( output_file.fail() ) throw runtime_error("Could not open filename '" + output + "'");
    factory.output_feature_matrix(output_file, cache, sample);

    // Print feature metadata
    output = options.workspace_ + "/feature-info.io";
    ofstream infofile(output);
    if( infofile.fail() ) throw runtime_error("Could not open filename '" + output + "'");
    factory.output_feature_info(infofile, cache, sample);
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

    const SLTP::DL::Sample sample = parse_input_sample(options.workspace_ + "/sample.io");
    const std::vector<std::string> nominals = parse_nominals(options.workspace_ + "/nominals.io");
#if 0
    cout << "SAMPLE: #objects=" << sample->num_objects()
         << ", #predicates=" << sample->num_predicates()
         << ", #grounded-predicates=" << sample->num_grounded_predicates()
         << ", #states=" << sample->num_states()
         << endl;
#endif
    SLTP::DL::Factory factory("test", nominals, options.complexity_bound_);

    SLTP::DL::Cache cache;
    factory.generate_basis(sample);
    factory.generate_roles(cache, sample);
    factory.generate_concepts(cache, sample);
    factory.generate_features(cache, sample);
//    factory.report_dl_data(cout);
    factory.log_all_concepts_and_features(cache, sample, options.workspace_);

    output_results(options, factory, sample, cache);

    return 0;
}

