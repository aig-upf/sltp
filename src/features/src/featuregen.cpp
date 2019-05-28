
#include <iostream>
#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>

#include <sltp/features.hxx>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
namespace po = boost::program_options;


//! Command-line option processing
struct Options {
    std::string workspace;
    unsigned complexity_bound;
    unsigned dist_complexity_bound;
    unsigned cond_complexity_bound;
    bool print_denotations;

    Options(int argc, const char **argv) {
        po::options_description description("Generate a pool of non-redundante candidate features given"
                                            " a sample of states.\n\nOptions");

        description.add_options()
                ("help,h", "Display this help message and exit.")
                ("verbose,v", "Show extra debugging messages.")

                ("workspace,w", po::value<std::string>()->required(),
                 "Directory where the input and output files will be expected / left.")

                ("complexity-bound", po::value<unsigned>()->required(),
                        "Maximum feature complexity for standard features.")
                ("dist-complexity-bound", po::value<unsigned>()->required(),
                        "Maximum feature complexity for distance features.")
                 ("cond-complexity-bound", po::value<unsigned>()->required(),
                         "Maximum feature complexity for conditional features.")

                ("print-denotations", "Whether print the denotations of all generated features.")
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
        complexity_bound = vm["complexity-bound"].as<unsigned>();
        dist_complexity_bound = vm["dist-complexity-bound"].as<unsigned>();
        cond_complexity_bound = vm["cond-complexity-bound"].as<unsigned>();
        print_denotations = vm.count("print-denotations") > 0;
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
    string output(options.workspace + "/feature-matrix.io");
    ofstream output_file(output);
    if( output_file.fail() ) throw runtime_error("Could not open filename '" + output + "'");
    factory.output_feature_matrix(output_file, cache, sample);

    // Print feature metadata
    output = options.workspace + "/feature-info.io";
    ofstream infofile(output);
    if( infofile.fail() ) throw runtime_error("Could not open filename '" + output + "'");
    factory.output_feature_info(infofile, cache, sample);
}

int main(int argc, const char **argv) {
    Options options(argc, argv);

    const SLTP::DL::Sample sample = parse_input_sample(options.workspace + "/sample.io");
    std::vector<std::string> nominals = parse_nominals(options.workspace + "/nominals.io");
#if 0
    cout << "SAMPLE: #objects=" << sample->num_objects()
         << ", #predicates=" << sample->num_predicates()
         << ", #grounded-predicates=" << sample->num_grounded_predicates()
         << ", #states=" << sample->num_states()
         << endl;
#endif
    SLTP::DL::Factory factory(nominals, options.complexity_bound,
            options.dist_complexity_bound, options.cond_complexity_bound);

    SLTP::DL::Cache cache;
    factory.generate_basis(sample);
    // TODO Make this optional adding an extra command-line option:
    auto forced_goal_features = factory.generate_goal_concepts_and_roles(cache, sample);

    factory.generate_roles(cache, sample);

    auto concepts = factory.generate_concepts(cache, sample);
    factory.generate_features(concepts, cache, sample, forced_goal_features);
//    factory.report_dl_data(cout);
    factory.log_all_concepts_and_features(concepts, cache, sample, options.workspace, options.print_denotations);

    output_results(options, factory, sample, cache);

    return 0;
}

