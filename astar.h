//
// Created by egorb on 18.01.2023.
//

#ifndef KP_DA_ASTAR_H
#define KP_DA_ASTAR_H

#include "utils.h"
#include "queue"

std::vector<node> readNodes(std::ifstream &graph) {
    uint numberOfNodes;
    graph >> numberOfNodes;
    std::vector<node> V(numberOfNodes);
    for (uint i = 0; i < numberOfNodes; ++i) {
        graph >> V[i];
    }

    return V;
}

uint getPos(uint id, const std::vector<node> &V) {
    return std::lower_bound(V.begin(), V.end(), id) - V.begin();
}

std::vector<uint> getNeighbors(uint pos, std::ifstream &graph, const uint offset, const std::vector<node> &V) {
    graph.seekg(offset + V[pos].offset);
    int end = (pos + 1 < V.size() ? (int) V[pos + 1].offset : -1);
    std::vector<uint> neighbors;
    uint to;
    while (graph.tellg() != offset + end && graph >> to) {
        neighbors.push_back(getPos(to, V));
    }

    return neighbors;
}

void find(std::ifstream &graph, const uint offset, const std::vector<node> &V, const uint start, const uint finish,
          std::ofstream &output) {
    if (start == finish) {
        output << "0\n";
        return;
    }
    uint finishPos = getPos(finish, V);

    auto h = [&V, finishPos](uint pos) -> double {
        static node finishNode = V[finishPos];

        return dist(V[pos], finishNode);
    };

    std::vector<double> d(V.size(), -1);
    std::priority_queue<std::pair<double, uint>, std::vector<std::pair<double, uint>>, std::greater<>> pq;
    uint startPos = getPos(start, V);
    d[startPos] = 0;
    pq.emplace(d[startPos] + h(startPos), startPos);
    while (!pq.empty()) {

        auto [w, source] = pq.top();
        pq.pop();
        if (source == finishPos) {
            output << d[source] << "\n";
            return;
        }
        if (w > d[source] + h(source) + eps) {
            continue;
        }
        std::vector<uint> neighbors = getNeighbors(source, graph, offset, V);

        for (const auto &destination: neighbors) {
            double weight = dist(V[source], V[destination]);
            if (fabs(d[destination] + 1) < eps || d[destination] > d[source] + weight + eps) {
                d[destination] = d[source] + weight;
                pq.push({d[destination] + h(destination), destination});
            }
        }
    }
    output << "-1\n";
}

void solve(const std::string &graphFilename, const std::string &inputFilename, const std::string &outputFilename,
           bool isDetailed = false) {
    std::ifstream graph(graphFilename);
    std::ifstream input(inputFilename);
    std::ofstream output(outputFilename);

    std::vector<node> V = readNodes(graph);
    uint offset = graph.tellg();

    uint start, finish;
    while (input >> start >> finish) {
        find(graph, offset, V, start, finish, output);
    }
}

#endif //KP_DA_ASTAR_H
