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
#include <cnf/transition_classification.h>
#include <cnf/aaai19_generator.h>
#include <common/helpers.h>


using namespace std;
namespace po = boost::program_options;

//! Command-line option processing
sltp::cnf::Options parse_options(int argc, const char **argv) {
    sltp::cnf::Options options;

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

        ("enforce-features,e", po::value<std::string>()->default_value(""),
         "Comma-separated (no spaces!) list of IDs of feature we want to enforce in the abstraction.")


        ("v_slack", po::value<double>()->default_value(2),
         "The slack value for the maximum allowed value for V_pi(s) = slack * V^*(s)")

        ("use-incremental-refinement",
         "In the transition-separation encoding, whether to use the incremental refinement approach")

        ("distinguish-goals",
         "In the transition-separation encoding, whether to post constraints to ensure distinguishability of goals")

        ("encoding", po::value<std::string>()->default_value("basic"),
             "The encoding to be used (options: {basic, d2tree, separation}).")

        ("use-equivalence-classes",
         "In the transition-separation encoding, whether we want to exploit the equivalence relation "
         "among transitions given by the feature pool")



        ("use-feature-dominance",
         "In the transition-separation encoding, whether we want to exploit the dominance among features to ignore "
         "dominated features and reduce the size of the encoding.")

        ("distinguish-transitions-locally",
         "In the transition-separation CNF encoding, whether to distinguish good transitions *only from*"
         " unmarked transitions starting in the same state")
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

    options.workspace = vm["workspace"].as<std::string>();
    options.prune_redundant_states = vm.count("prune-redundant-states") > 0;
    options.verbose = vm.count("verbose") > 0;
    options.distinguish_transitions_locally = vm.count("distinguish-transitions-locally") > 0;
    options.use_equivalence_classes = vm.count("use-equivalence-classes") > 0;
    options.use_feature_dominance = vm.count("use-feature-dominance") > 0;
    options.use_incremental_refinement = vm.count("use-incremental-refinement") > 0;
    options.distinguish_goals = vm.count("distinguish-goals") > 0;
    options.v_slack = vm["v_slack"].as<double>();

    auto enc = vm["encoding"].as<std::string>();
    if (enc == "basic") options.encoding = sltp::cnf::Options::Encoding::Basic;
    else if (enc == "d2tree") options.encoding = sltp::cnf::Options::Encoding::D2Tree;
    else if  (enc == "separation") options.encoding = sltp::cnf::Options::Encoding::TransitionSeparation;
    else throw po::validation_error(po::validation_error::invalid_option_value, "encoding");


    // Split the comma-separated list of enforced feature IDS
    auto enforced_str = vm["enforce-features"].as<std::string>();
//        std::cout << "\nenforced_str: " << enforced_str << std::endl;
    if (!enforced_str.empty()) {
        std::stringstream ss(enforced_str);
        while (ss.good()) {
            std::string substr;
            getline(ss, substr, ',');
            if (!substr.empty()) {
                options.enforced_features.push_back((unsigned) atoi(substr.c_str()));
            }
        }
    }
//        for (auto x:enforced_features) std::cout << "\n" << x << std::endl;

    return options;
}

sltp::cnf::CNFGenerationOutput
write_encoding(CNFWriter& wr, const sltp::TrainingSet& sample, const sltp::cnf::Options& options) {
    if (options.use_separation_encoding()) {
        sltp::cnf::TransitionClassificationEncoding generator(sample, options);
        return generator.refine_theory(wr);

    } else {
        if (options.prune_redundant_states) {
            // If indicated by the user, prune those states that appear redundant for the given feature pool
            auto resample = sltp::cnf::AAAI19Generator::preprocess_sample(sample, options);
            sltp::cnf::AAAI19Generator gen(resample, options);
            return gen.write(wr);
        } else {
            sltp::cnf::AAAI19Generator gen(sample, options);
            return gen.write(wr);
        }
    }
}

int main(int argc, const char **argv) {
    float start_time = Utils::read_time_in_seconds();
    auto options = parse_options(argc, argv);

    // Read input training set
    sltp::TrainingSet trset(
            read_feature_matrix(options.workspace, options.verbose),
            read_transition_data(options.workspace, options.verbose),
            read_input_sample(options.workspace));

    std::cout << "Training sample: " << trset << std::endl;

    // We write the MaxSAT instance progressively as we generate the CNF. We do so into a temporary "*.tmp" file
    // which will be later processed by the Python pipeline to inject the value of the TOP weight, which we can
    // know only when we finish writing all clauses

    auto wsatstream = get_ofstream(options.workspace + "/theory.wsat.tmp");
    auto allvarsstream = get_ofstream(options.workspace + "/allvars.wsat");

    CNFWriter writer(wsatstream, &allvarsstream);
    auto output = write_encoding(writer, trset, options);

    wsatstream.close();
    allvarsstream.close();

    if (output != sltp::cnf::CNFGenerationOutput::ValidationCorrectNoRefinementNecessary) {
        // Write some characteristics of the CNF to a different file
        auto topstream = get_ofstream(options.workspace + "/top.dat");
        topstream << writer.top() << " " << writer.nvars() << " " << writer.nclauses();
        topstream.close();

        float total_time = Utils::read_time_in_seconds() - start_time;
        std::cout << "Total-time: " << total_time << std::endl;
        std::cout << "CNF Theory: " << writer.nvars() << " vars + " << writer.nclauses() << " clauses" << std::endl;
    }

    if (output == sltp::cnf::CNFGenerationOutput::UnsatTheory) {
        std::cout << Utils::warning() << "CNF theory is UNSAT" << std::endl;
    }

    return static_cast<std::underlying_type_t<sltp::cnf::CNFGenerationOutput>>(output);
}
