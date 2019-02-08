
#include <sltp/features.hxx>
#include <sltp/io.hxx>
#include <iostream>


struct Options {
    std::string worskspace;
};

Options process_input(int argc, const char **argv);

void print_usage(std::ostream &os, const std::string &execname);


int main(int argc, const char **argv) {
    auto options = process_input(argc, argv);

    // Read primitive concepts and files from input file
    auto primitive_denotations = sltp::io::read_denotation_matrix(options.worskspace + "/primitives.csv");
    auto static_denotations = sltp::io::read_denotation_matrix(options.worskspace + "/statics.csv");


    // Generate concepts based on grammar, pruning duplicate concepts


    // Generate Features and prune duplicates


    // Generate output "state x feature" matrix
    // We will likely need to output somehow the structure of the (non-redundant) concepts and features,
    // so that they can be reconstructed from python.




    return 0;
}



Options process_input(int argc, const char **argv) {

    // Print call
    std::cout << "Call:";
    for(unsigned i = 0; i < argc; ++i) std::cout << " " << argv[i];
    std::cout << std::endl;

    std::string execname(*argv); // executable name
    if(argc != 3) { // check we have enough arguments
        print_usage(std::cout, execname);
        exit(1);
    }

    // read arguments
    Options options;
    options.worskspace = std::string(argv[1]);
    std::cout << "Using workspace directory " << options.worskspace << std::endl;
    return options;
}

void print_usage(std::ostream& os, const std::string& execname) {
    std::cout << std::endl
         << "Usage: " << execname << " <input-file> <output-file>" << std::endl
         << std::endl
         << "where" << std::endl
         << "    <input-file> is the path to a CSV file containing the denotations of primitive concepts and roles "
            "                  from which to start the feature generation process." << std::endl
         << "    <output-file> is the path to the output CSV file where the denotation matrix of all generated features "
            "                  will be left." << std::endl
         << std::endl
         ;
    }
