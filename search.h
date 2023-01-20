//
// Created by egorb on 18.01.2023.
//

#ifndef KP_DA_ASTAR_CPP
#define KP_DA_ASTAR_CPP

#include <unordered_set>
#include "utils.h"
#include "queue"
#include "unordered_map"
#include "iomanip"


std::vector<compactNode<long long>> readNodes(std::ifstream &graph) {
    uint numberOfNodes;
    graph.read(reinterpret_cast<char *>(&numberOfNodes), sizeof(numberOfNodes));
    std::vector<compactNode<long long>> V(numberOfNodes);
    for (uint i = 0; i < numberOfNodes; ++i) {
        long long offset = graph.tellg();
        node v;
        graph.read(reinterpret_cast<char *>(&v.id), sizeof(v.id));
        graph.read(reinterpret_cast<char *>(&v.phi), sizeof(v.phi));
        graph.read(reinterpret_cast<char *>(&v.lambda), sizeof(v.lambda));
        graph.read(reinterpret_cast<char *>(&v.offset), sizeof(v.offset));
        V[i] = compactNode(v.id, offset);
    }

    return V;
}


uint getPos(uint id, const std::vector<compactNode<long long>> &V) {
    compactNode<long long> s(id, 0);
    return std::lower_bound(V.begin(), V.end(), s) - V.begin();
}

/*
 * Getting node by position
 */
node getNode(uint pos, std::ifstream &graph, const std::vector<compactNode<long long>> &V) {
    if (graph.tellg() == -1)graph.clear();
    graph.seekg(V[pos].offset, std::ifstream::beg);

    uint id, offset;
    double phi, lambda;
    graph.read(reinterpret_cast<char *>(&id), sizeof(id));
    graph.read(reinterpret_cast<char *>(&phi), sizeof(phi));
    graph.read(reinterpret_cast<char *>(&lambda), sizeof(lambda));
    graph.read(reinterpret_cast<char *>(&offset), sizeof(offset));

    if (id != V[pos].id) {
        assert(id == V[pos].id);
    }

    return {id, phi, lambda, offset};
}

std::vector<uint>
getNeighbors(uint pos, std::ifstream &graph, const long long offset, const std::vector<compactNode<long long>> &V) {

    long long end = (pos + 1 < V.size() ? (long long) getNode(pos + 1, graph, V).offset : -1) + offset;
    long long newPosition = offset + (long long) getNode(pos, graph, V).offset;

    if (graph.tellg() == -1)graph.clear();
    graph.seekg(newPosition, std::ifstream::beg);

    std::vector<uint> neighbors;
    uint to;
    while (graph.tellg() != end && graph.read(reinterpret_cast<char *>(&to), sizeof(to))) {
        neighbors.push_back(to);
    }

    return neighbors;
}


double h(std::ifstream &graph, const std::vector<compactNode<long long>> &V, uint startPos, node finishNode) {
    return dist(getNode(startPos, graph, V), finishNode);
}

void
find(std::ifstream &graph, const long long offset, const std::vector<compactNode<long long>> &V, const uint start,
     const uint finish,
     std::ofstream &output, bool isDetailed) {
    uint finishPos = getPos(finish, V);
    node finishNode = getNode(finishPos, graph, V);

    std::vector<double> d(V.size(), -1);
    std::vector<uint> prev(V.size());

    std::priority_queue<std::pair<double, uint>, std::vector<std::pair<double, uint>>, std::greater<>> pq;

    uint startPos = getPos(start, V);
    d[startPos] = 0;
    prev[startPos] = 0;
    pq.emplace(d[startPos] + h(graph, V, startPos, finishNode), startPos);

    while (!pq.empty()) {
        auto [w, source] = pq.top();
        pq.pop();

        if (source == finishPos) {
            break;
        }

        if (w > d[source] + h(graph, V, source, finishNode) + eps) {
            continue;
        }

        std::vector<uint> neighbors = getNeighbors(source, graph, offset, V);

        node sourceNode = getNode(source, graph, V);

        for (const auto &destination: neighbors) {
            double weight = dist(sourceNode, getNode(destination, graph, V));
            if (fabs(d[destination] + 1) < eps || d[destination] > d[source] + weight + eps) {
                d[destination] = d[source] + weight;
                prev[destination] = source;
                pq.emplace(d[destination] + h(graph, V, destination, finishNode), destination);
            }
        }
    }
    output << std::fixed << std::setprecision(7);
    if (fabs(d[finishPos] + 1) < eps) {
        output << "-1\n";
        if (isDetailed) {
            output << "0" << "\n";
        }
    } else {
        output << d[finishPos] << "\n";
        if (isDetailed) {
            d.clear();
            d.shrink_to_fit();

            std::vector<uint> path;
            for (uint cur = finishPos;; cur = prev[cur]) {
                path.push_back(cur);
                if (cur == startPos)break;
            }
            std::reverse(path.begin(), path.end());
            output << path.size() << " ";
            for (auto &pos: path) {
                pos = V[pos].id;
                output << pos << " ";
            }
            output << "\n";
        }
    }
}

void solve(const std::string &graphFilename, const std::string &inputFilename, const std::string &outputFilename,
           bool isDetailed = false) {
    std::ifstream graph(graphFilename, std::ios::binary);
    std::ifstream input(inputFilename);
    std::ofstream output(outputFilename);


    std::vector<compactNode<long long>> V = readNodes(graph);
    sort(V.begin(), V.end());

    long long offset = graph.tellg();

    uint start, finish;
    while (input >> start >> finish) {
        find(graph, offset, V, start, finish, output, isDetailed);
    }
}

#endif //KP_DA_ASTAR_CPP
