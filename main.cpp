#include <cstring>
#include "iostream"
#include "utils.h"
#include "string"
#include "search.h"

int main(int argc, char **argv) {
    if (argc < 8) {
        std::cout << "Not enough parameters" << "\n";
        return 1;
    }
    if (strcmp(argv[1], "preprocess") == 0) {
        std::string nodesFilename, edgesFilename, outputFilename;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--nodes") == 0) {
                nodesFilename = argv[i + 1];
            } else if (strcmp(argv[i], "--edges") == 0) {
                edgesFilename = argv[i + 1];
            } else if (strcmp(argv[i], "--output") == 0) {
                outputFilename = argv[i + 1];
            }
        }
        if (nodesFilename.empty() || edgesFilename.empty() || outputFilename.empty()) {
            std::cout << "Wrong parameters" << "\n";
            return 1;
        }
        prepareFile(nodesFilename, edgesFilename, outputFilename);
    } else if (strcmp(argv[1], "search") == 0) {

        std::string graphFilename, inputFilename, outputFilename;
        bool isDetailed = false;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--graph") == 0) {
                graphFilename = argv[i + 1];
            } else if (strcmp(argv[i], "--input") == 0) {
                inputFilename = argv[i + 1];
            } else if (strcmp(argv[i], "--output") == 0) {
                outputFilename = argv[i + 1];
            } else if (strcmp(argv[i], "--full-output") == 0) {
                isDetailed = true;
            }
        }
        if (graphFilename.empty() || inputFilename.empty() || outputFilename.empty()) {
            std::cout << "Wrong parameters" << "\n";
            return 1;
        }

        solve(graphFilename, inputFilename, outputFilename, isDetailed);
    } else {
        std::cout << "Wrong mode\n";
        return 1;
    }
//    prepareFile("n_gen.txt", "e_gen.txt", "g.txt");
//    solve("g.txt", "q_gen.txt", "res.txt", true);
    return 0;
}