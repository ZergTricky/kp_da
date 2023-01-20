//
// Created by egorb on 18.01.2023.
//
#ifndef KP_DA_UTILS_H
#define KP_DA_UTILS_H

#include <chrono>
#include "iostream"
#include "fstream"
#include "cmath"
#include "cassert"
#include "vector"
#include "string"


using uint = unsigned int;

const double PI = acos(-1);
const double eps = 1e-9;
const double R = 6371;

inline double rad(double angle) {
    return PI / 180.0 * angle;
}

class node {
public:
    uint id{};
    uint offset = UINT32_MAX;

    node() = default;

    node(uint id) : id(id) {};

    node(uint id, double phi, double lambda) : node(id, phi, lambda, 0) {}

    node(uint id, double phi, double lambda, uint offset) : id(id), offset(offset), phi(rad(phi)),
                                                            lambda(rad(lambda)) {}

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

    double phi = 0;
    double lambda = 0;

    friend double dist(const node &, const node &);
};

template<typename T>
class compactNode {
public:
    uint id = 0;
    T offset = UINT32_MAX;

    compactNode() = default;

    compactNode(uint id) : id(id) {}

    compactNode(uint id, T offset) : id(id), offset(offset) {}

    bool operator<(const compactNode &other) const {
        if (id != other.id) {
            return id < other.id;
        }
        return offset < other.offset;
    }
};

double dist(const node &a, const node &b) {
    double cosd = sin(a.phi) * sin(b.phi) + cos(a.phi) * cos(b.phi) * cos(a.lambda - b.lambda);
    assert(fabs(cosd) < 1 + eps);
    double d = acos(cosd);
    return R * d;
}

using edge = std::pair<uint, uint>;

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


uint countNodes(std::ifstream &nodes) {
    uint cnt = 0;
    uint id;
    double phi, lambda;
    while (nodes >> id >> phi >> lambda) {
        ++cnt;
    }

    nodes.clear();
    nodes.seekg(0);

    return cnt;
}

std::vector<compactNode<uint>> proceedNodes(std::ifstream &nodes) {
    std::vector<compactNode<uint>> V;

    uint id;
    double phi, lambda;
    while (nodes >> id >> phi >> lambda) {
        V.emplace_back(id);
    }
    V.shrink_to_fit();

    return V;
}

uint countEdges(std::ifstream &edges) {
    uint cnt = 0;
    uint k;
    while (edges >> k) {
        for (uint i = 0; i < k; ++i) {
            uint _;
            edges >> _;
            if (i > 0) {
                cnt += 2;
            }
        }
    }
    edges.clear();
    edges.seekg(0);

    return cnt;
}

std::vector<edge> proceedEdges(std::ifstream &edges) {
    std::vector<edge> E;

    uint k;
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

    return E;
}

uint getNodePos(uint id, const std::vector<compactNode<uint>> &V) {
    compactNode<uint> s(id, 0);
    return std::lower_bound(V.begin(), V.end(), s) - V.begin();
}

void setOffset(std::ofstream &output, std::vector<compactNode<uint>> &V, const std::vector<edge> &E) {
    if (E.empty()) {
        return;
    }
    uint cur = 0;
    uint i;
    for (i = 0; i < E.size(); ++i) {
        if (E[i].first != E[cur].first) {
            uint pos = getNodePos(E[cur].first, V);
            V[pos].offset = output.tellp();
            for (uint j = cur; j < i; ++j) {
                uint p = getNodePos(E[j].second, V);
                output.write(reinterpret_cast<const char *>(&p), sizeof(p));
            }
            cur = i;
        }
    }
    uint pos = getNodePos(E[cur].first, V);
    V[pos].offset = output.tellp();
    for (uint j = cur; j < i; ++j) {
        uint p = getNodePos(E[j].second, V);
        output.write(reinterpret_cast<const char *>(&p), sizeof(p));
    }
}

void printEdges(std::ofstream &output, const std::vector<edge> &E, const std::vector<compactNode<uint>> &V) {
    if (E.empty()) {
        return;
    }
    uint cur = 0;
    uint i;
    for (i = 0; i < E.size(); ++i) {
        if (E[i].first != E[cur].first) {
            for (uint j = cur; j < i; ++j) {
                auto p = getNodePos(E[j].second, V);
                output.write(reinterpret_cast<const char *>(&p), sizeof(p));
            }
            cur = i;
        }
    }
    for (uint j = cur; j < i; ++j) {
        uint p = getNodePos(E[j].second, V);
        output.write(reinterpret_cast<const char *>(&p), sizeof(p));
    }
}

void printNodes(std::ofstream &output, std::ifstream &nodes, const std::vector<compactNode<uint>> &V) {
    if (nodes.tellg() == -1)nodes.clear();
    nodes.seekg(0);
    uint id;
    double phi, lambda;
    while (nodes >> id >> phi >> lambda) {
        compactNode<uint> s(id, 0);
        auto iter = std::lower_bound(V.begin(), V.end(), s);
        assert(iter->id == id);
        output.write(reinterpret_cast<const char *>(&iter->id), sizeof(iter->id));
        output.write(reinterpret_cast<const char *>(&phi), sizeof(phi));
        output.write(reinterpret_cast<const char *>(&lambda), sizeof(lambda));
        output.write(reinterpret_cast<const char *>(&iter->offset), sizeof(iter->offset));
    }
}

void
prepareFile(const std::string &nodesFilename, const std::string &edgesFilename, const std::string &outputFilename) {
    std::ifstream nodes(nodesFilename, std::ios::in);
    std::ifstream edges(edgesFilename, std::ios::in);
    std::ofstream output(outputFilename, std::ios::binary | std::ios::out);

    std::vector<compactNode<uint>> V = proceedNodes(nodes);
    std::sort(V.begin(), V.end());
    std::vector<edge> E = proceedEdges(edges);
    std::sort(E.begin(), E.end(), [](const edge &lhs, const edge &rhs) -> bool {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });

    setOffset(output, V, E);

    for (int i = (int) V.size() - 2; i >= 0; --i) {
        V[i].offset = std::min(V[i].offset, V[i + 1].offset);
    }

    // Clear output file for correct offset
    output.close();
    output.open(outputFilename, std::ios::binary | std::ios::out);
    uint size = V.size();
    output.write(reinterpret_cast<const char *>(&size), sizeof(size));

    printNodes(output, nodes, V);
    printEdges(output, E, V);
}

#endif //KP_DA_UTILS_H