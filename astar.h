//
// Created by egorb on 18.01.2023.
//

#ifndef KP_DA_ASTAR_H
#define KP_DA_ASTAR_H

#include "utils.h"
#include "queue"
#include "unordered_map"
#include "unordered_set"


std::vector<compactNode> readNodes(std::ifstream &graph) {
    uint numberOfNodes;
    graph >> numberOfNodes;
    std::vector<compactNode> V(numberOfNodes);
    for (uint i = 0; i < numberOfNodes; ++i) {
        long long offset = graph.tellg();
        node v;
        graph >> v;
        V[i] = compactNode(v.id, offset); // todo won't fit
    }

    return V;
}

uint getPos(uint id, const std::vector<compactNode> &V) {
    return std::lower_bound(V.begin(), V.end(), id) - V.begin();
}

/*
 * Getting node by position
 */
node getNode(uint pos, std::ifstream &graph, const std::vector<compactNode> &V) {
    if (graph.tellg() == -1)graph.clear();
    graph.seekg(V[pos].offset);
    uint id, offset;
    double phi, lambda;
    graph >> id >> phi >> lambda >> offset;
    assert(id == V[pos].id);
    return {id, phi, lambda, offset};
}

std::vector<uint>
getNeighbors(uint pos, std::ifstream &graph, const long long offset, const std::vector<compactNode> &V) {

    long long end = (pos + 1 < V.size() ? (long long) getNode(pos + 1, graph, V).offset : -1) + offset;
    long long newPosition = (long long) offset + (long long) getNode(pos, graph, V).offset;

    if (graph.tellg() == -1)graph.clear();
    graph.seekg(newPosition);
    assert(graph.tellg() == newPosition);

    std::vector<uint> neighbors;
    uint to;
    while (graph.tellg() != end && graph >> to) {
        neighbors.push_back(getPos(to, V));
    }

    return neighbors;
}

struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

void
find(std::ifstream &graph, const long long offset, const std::vector<compactNode> &V, const uint start,
     const uint finish,
     std::ofstream &output) {
#ifdef DEBUG
    std::cout << start << " " << finish << std::endl;
#endif
    if (start == finish) {
        output << "0\n";
        return;
    }
    uint finishPos = getPos(finish, V);

    auto h = [&graph, &V, finishPos](uint pos) -> double {
        static node finishNode = getNode(finishPos, graph, V);

        return dist(getNode(pos, graph, V), finishNode);
    };

    // todo Better use unordered_map?
    std::unordered_map<uint, double, custom_hash> d;
    std::unordered_set<uint, custom_hash> used;
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
        if (used.contains(source)) {
            continue;
        }
        used.insert(source);
//        if (used.size() == 200000)exit(0);
        std::vector<uint> neighbors = getNeighbors(source, graph, offset, V);

        node sourceNode = getNode(source, graph, V);

        for (const auto &destination: neighbors) {
            double weight = dist(sourceNode, getNode(destination, graph, V));
            if (!d.contains(destination) || d[destination] > d[source] + weight + eps) {
                d[destination] = d[source] + weight;
                pq.emplace(d[destination] + h(destination), destination);
            }
        }
    }
    output << "-1\n";
}

void solve(const std::string &graphFilename, const std::string &inputFilename, const std::string &outputFilename,
           bool isDetailed = false) {
    std::ifstream graph(graphFilename, std::ios::binary);
    std::ifstream input(inputFilename);
    std::ofstream output(outputFilename, std::ios::binary);

#ifdef DEBUG
    auto startNodes = std::chrono::steady_clock::now();
#endif
    std::vector<compactNode> V = readNodes(graph);
    sort(V.begin(), V.end());
#ifdef DEBUG
    std::cout << "Nodes done! in " << since(startNodes).count() / 1000 << "sec " << std::endl;
#endif
    long long offset = graph.tellg();

    uint start, finish;
    while (input >> start >> finish) {
        find(graph, offset, V, start, finish, output);
    }
}

#endif //KP_DA_ASTAR_H
