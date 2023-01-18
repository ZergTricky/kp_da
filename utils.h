//
// Created by egorb on 18.01.2023.
//

#ifndef KP_DA_UTILS_H
#define KP_DA_UTILS_H

#include <chrono>
#include "fstream"
#include "cmath"
#include "cassert"
#include "vector"

#include "iostream"

#define DEBUG

using uint = unsigned int;

const double PI = acos(-1);
const double eps = 1e-9;
const double R = 6371;

inline double rad(double angle) {
    return PI / 180.0 * angle;
}

class node {
public:
    uint id;
    uint offset = 0;

    node() = default;

    node(uint id) : id(id) {};

    node(uint id, double phi, double lambda) : id(id), phi(rad(phi)), lambda(rad(lambda)) {}

    ~node() = default;

    bool operator<(const node &other) const {
        return id < other.id;
    }

    friend std::istream &operator>>(std::istream &is, node &nd) {
        is >> nd.id >> nd.phi >> nd.lambda >> nd.offset;
        return is;
    }

    friend std::ostream &operator<<(std::ostream &os, const node &nd) {
        os << nd.id << " " << nd.phi << " " << nd.lambda << " " << nd.offset << "\n";
        return os;
    }

private:
    double phi = 0;
    double lambda = 0;

    friend double dist(node, node);
};

using edge = std::pair<uint, uint>;

double dist(node a, node b) {
    double cosd = sin(a.phi) * sin(b.phi) + cos(a.phi) * cos(b.phi) * cos(a.lambda - b.lambda);
    assert(fabs(cosd) < 1 + eps);
    double d = acos(cosd);
    return R * d;
}

std::vector<node> proceedNodes(std::ifstream &nodes) {
    std::vector<node> V;
    uint id;
    double phi, lambda;
    while (nodes >> id >> phi >> lambda) {
        V.emplace_back(id, phi, lambda);
    }
    V.shrink_to_fit();

    return V;
}

std::vector<edge> proceedEdges(std::ifstream &edges) {
    std::vector<edge> E;
    int k;
    while (edges >> k) {
        uint last;
        for (uint i = 0; i < k; ++i) {
            uint cur;
            edges >> cur;
            if (i > 0) {
                E.emplace_back(last, cur);
                E.emplace_back(cur, last);
            }
            last = cur;
        }
    }
    E.shrink_to_fit();

    return E;
}

void setOffset(std::ofstream &output, std::vector<node> &V, const std::vector<edge> &E) {
    uint last;
    std::vector<uint> cur;
    for (uint i = 0; i < E.size(); ++i) {
        if (i == 0 || E[i].first != last) {
            if (!cur.empty()) {
                auto iter = std::lower_bound(V.begin(), V.end(), last);
                iter->offset = output.tellp();

                for (uint j = 0; j < cur.size(); ++j) {
                    output << cur[j];
                    if (j + 1 < cur.size())output << " ";
                }
                output << "\n";
            }
            cur.clear();
        }
        cur.emplace_back(E[i].second);

        last = E[i].first;
    }
    if (!cur.empty()) {
        auto iter = std::lower_bound(V.begin(), V.end(), last);
        iter->offset = output.tellp();

        for (uint j = 0; j < cur.size(); ++j) {
            output << cur[j];
            if (j + 1 < cur.size())output << " ";
        }
        output << "\n";
    }
}

void printEdges(std::ofstream &output, const std::vector<edge> &E) {
    uint last;
    std::vector<uint> cur;
    for (uint i = 0; i < E.size(); ++i) {
        if (i == 0 || E[i].first != last) {
            if (!cur.empty()) {
                for (uint j = 0; j < cur.size(); ++j) {
                    output << cur[j];
                    if (j + 1 < cur.size())output << " ";
                }
                output << "\n";
            }
            cur.clear();
        }
        cur.emplace_back(E[i].second);

        last = E[i].first;
    }
    if (!cur.empty()) {
        for (uint j = 0; j < cur.size(); ++j) {
            output << cur[j];
            if (j + 1 < cur.size())output << " ";
        }
        output << "\n";
    }
}

void printNodes(std::ofstream &output, const std::vector<node> &V) {
    for (const auto &nd: V) {
        output << nd;
    }
}

#ifdef DEBUG

template<
        class result_t   = std::chrono::milliseconds,
        class clock_t    = std::chrono::steady_clock,
        class duration_t = std::chrono::milliseconds
>
auto since(std::chrono::time_point<clock_t, duration_t> const &start) {
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

#endif


void
prepareFile(const std::string &nodesFilename, const std::string &edgesFilename, const std::string &outputFilename) {


    std::ifstream nodes;
    nodes.open(nodesFilename);

    std::ifstream edges;
    edges.open(edgesFilename);

    std::ofstream output;
    output.open(outputFilename, std::ios::binary);

#ifdef DEBUG
    auto startNodes = std::chrono::steady_clock::now();
#endif
    std::vector<node> V = proceedNodes(nodes);
    std::sort(V.begin(), V.end());
#ifdef DEBUG
    std::cout << "Nodes done! in " << since(startNodes).count() / 1000 << "sec " << std::endl;
#endif


#ifdef DEBUG
    auto startEdges = std::chrono::steady_clock::now();
#endif
    std::vector<edge> E = proceedEdges(edges);

    std::sort(E.begin(), E.end(), [](const edge &lhs, const edge &rhs) -> bool {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });
#ifdef DEBUG
    std::cout << "Edges done! in " << since(startEdges).count() / 1000 << "sec " << std::endl;
#endif
    setOffset(output, V, E);

    // Clear output file for correct offset
    output.close();
    output.open(outputFilename);

    output << V.size() << "\n";
    printNodes(output, V);
    printEdges(output, E);
}

#endif //KP_DA_UTILS_H