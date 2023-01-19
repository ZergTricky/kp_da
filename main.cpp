#include "iostream"
#include "utils.h"
#include "astar.h"

int main(int argc, char **argv) {
//    prepareFile("n.txt", "e.txt", "g.txt");
    solve("g_europe.txt", "query.txt", "res.txt");
    return 0;
}