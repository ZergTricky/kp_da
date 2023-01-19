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


std::vector<compactNode<long long>> readNodes(FILE *graph) {
    uint numberOfNodes;
    fscanf(graph, "%u", &numberOfNodes);
    std::vector<compactNode<long long>> V(numberOfNodes);
    for (uint i = 0; i < numberOfNodes; ++i) {
        long long offset = ftell(graph);
        node v;
        fscanf(graph, "%u %lf %lf %u", &v.id, &v.phi, &v.lambda, &v.offset);
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
node getNode(uint pos, FILE *graph, const std::vector<compactNode<long long>> &V) {

    fseek(graph, V[pos].offset, 0);

    uint id, offset;
    double phi, lambda;
    fscanf(graph, "%u%lf%lf%u", &id, &phi, &lambda, &offset);

    if (id != V[pos].id) {
        assert(id == V[pos].id);
    }

    return {id, phi, lambda, offset};
}

std::vector<uint>
getNeighbors(uint pos, FILE *graph, const long long offset, const std::vector<compactNode<long long>> &V) {

    long long end = (pos + 1 < V.size() ? (long long) getNode(pos + 1, graph, V).offset : -1) + offset;
    long long newPosition = offset + (long long) getNode(pos, graph, V).offset;

    fseek(graph, newPosition, 0);

    std::vector<uint> neighbors;
    uint to;
    while (ftell(graph) != end && fscanf(graph, "%u", &to) > 0) {
        neighbors.push_back(to);
    }

    return neighbors;
}


double h(FILE *graph, const std::vector<compactNode<long long>> &V, uint startPos, node finishNode) {
    return dist(getNode(startPos, graph, V), finishNode);
}

void
find(FILE *graph, const long long offset, const std::vector<compactNode<long long>> &V, const uint start,
     const uint finish,
     FILE *output, bool isDetailed) {

#ifdef DEBUG
    std::cout << "Searching: " << start << " " << finish << std::endl;
#endif
    uint finishPos = getPos(finish, V);
    node finishNode = getNode(finishPos, graph, V);

    std::unordered_map<uint, double, custom_hash> d;
    std::unordered_map<uint, uint, custom_hash> prev;
    std::unordered_set<uint> s;

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
        s.insert(source);

        std::vector<uint> neighbors = getNeighbors(source, graph, offset, V);
        assert(neighbors.size() < 50);

        node sourceNode = getNode(source, graph, V);

        for (const auto &destination: neighbors) {
            double weight = dist(sourceNode, getNode(destination, graph, V));
            if (!d.contains(destination) || d[destination] > d[source] + weight + eps) {
                d[destination] = d[source] + weight;
                prev[destination] = source;
                pq.emplace(d[destination] + h(graph, V, destination, finishNode), destination);
            }
        }
    }
    if (!d.contains(finishPos)) {
        fprintf(output, "-1\n");
        if (isDetailed) {
            fprintf(output, "0\n");
        }
    } else {
        fprintf(output, "%.15f\n", d[finishPos]);
        if (isDetailed) {
            std::vector<uint> path;
            for (uint cur = finishPos;; cur = prev[cur]) {
                path.push_back(cur);
                if (cur == startPos)break;
            }
            std::reverse(path.begin(), path.end());
            fprintf(output, "%u ", (uint) path.size());
            for (auto &pos: path) {
                pos = V[pos].id;
                fprintf(output, "%u ", pos);
            }
            fprintf(output, "\n");
        }
    }
}

void solve(const std::string &graphFilename, const std::string &inputFilename, const std::string &outputFilename,
           bool isDetailed = false) {
    FILE *graph = fopen(graphFilename.c_str(), "rb");
    FILE *input = fopen(inputFilename.c_str(), "r");
    FILE *output = fopen(outputFilename.c_str(), "wb");

#ifdef DEBUG
    auto startNodes = std::chrono::steady_clock::now();
#endif

    std::vector<compactNode<long long>> V = readNodes(graph);
    sort(V.begin(), V.end());

#ifdef DEBUG
    std::cout << "Nodes done in " << since(startNodes).count() / 1000 << " sec" << std::endl;
#endif

    long long offset = ftell(graph);

    uint start, finish;
    while (fscanf(input, "%u %u", &start, &finish) > 0) {
        find(graph, offset, V, start, finish, output, isDetailed);
    }

    fclose(graph);
    fclose(input);
    fclose(output);
}

#endif //KP_DA_ASTAR_CPP
