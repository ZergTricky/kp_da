#include "iostream"
#include "utils.h"
#include "search.h"

int main(int argc, char **argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    prepareFile("n_gen.txt", "e_gen.txt", "g.txt");
    solve("g.txt", "q_gen.txt", "res.txt", true);
    return 0;
}