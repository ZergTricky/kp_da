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


uint countNodes(FILE *nodes) {
    uint cnt = 0;
    uint id;
    double phi, lambda;
    while (fscanf(nodes, "%u%lf%lf", &id, &phi, &lambda) > 0) {
        ++cnt;
    }

    fseek(nodes, 0, 0);

    return cnt;
}

std::vector<compactNode<uint>> proceedNodes(FILE *nodes) {
    uint N = countNodes(nodes);
    uint next = 0;
    std::vector<compactNode<uint>> V(N);

    uint id;
    double phi, lambda;
    while (fscanf(nodes, "%u%lf%lf", &id, &phi, &lambda) > 0) {
        V[next++] = id;
    }

    return V;
}

uint countEdges(FILE *edges) {
    uint cnt = 0;
    uint k;
    while (fscanf(edges, "%u", &k) > 0) {
        for (uint i = 0; i < k; ++i) {
            uint _;
            fscanf(edges, "%u", &_);
            if (i > 0) {
                cnt += 2;
            }
        }
    }
    fseek(edges, 0, 0);

    return cnt;
}

std::vector<edge> proceedEdges(FILE *edges) {
    uint N = countEdges(edges);
    uint next = 0;
    std::vector<edge> E(N);

    uint k;
    while (fscanf(edges, "%u", &k) > 0) {
        uint last;
        for (uint i = 0; i < k; ++i) {
            uint cur;
            fscanf(edges, "%u", &cur);
            if (i > 0) {
                E[next++] = {last, cur};
                E[next++] = {cur, last};
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

void setOffset(FILE *output, std::vector<compactNode<uint>> &V, const std::vector<edge> &E) {
    uint cur = 0;
    uint i;
    for (i = 0; i < E.size(); ++i) {
        if (E[i].first != E[cur].first) {
            uint pos = getNodePos(E[cur].first, V);
            V[pos].offset = ftell(output);
            for (uint j = cur; j < i; ++j) {
                uint p = getNodePos(E[j].second, V);
                fprintf(output, "%u", p);
                if (j + 1 < i) {
                    fprintf(output, " ");
                }
            }
            fprintf(output, "\n");
            cur = i;
        }
    }

    uint pos = getNodePos(E[cur].first, V);
    V[pos].offset = ftell(output);
    for (uint j = cur; j < i; ++j) {
        uint p = getNodePos(E[j].second, V);
        fprintf(output, "%u", p);
        if (j + 1 < i) {
            fprintf(output, " ");
        }
    }
    fprintf(output, "\n");
}

void printEdges(FILE *output, const std::vector<edge> &E, const std::vector<compactNode<uint>> &V) {
    uint cur = 0;
    uint i;
    for (i = 0; i < E.size(); ++i) {
        if (E[i].first != E[cur].first) {
            for (uint j = cur; j < i; ++j) {
                uint p = getNodePos(E[j].second, V);
                fprintf(output, "%u", p);
                if (j + 1 < i) {
                    fprintf(output, " ");
                }
            }
            fprintf(output, "\n");
            cur = i;
        }
    }
    for (uint j = cur; j < i; ++j) {
        uint p = getNodePos(E[j].second, V);
        fprintf(output, "%u", p);
        if (j + 1 < i) {
            fprintf(output, " ");
        }
    }
    fprintf(output, "\n");
}

void printNodes(FILE *output, FILE *nodes, const std::vector<compactNode<uint>> &V) {
    fseek(nodes, 0, 0);
    uint id;
    double phi, lambda;
    while (fscanf(nodes, "%u%lf%lf", &id, &phi, &lambda) > 0) {
        compactNode<uint> s(id, 0);
        auto iter = std::lower_bound(V.begin(), V.end(), s);
        assert(iter->id == id);
        fprintf(output, "%u %lf %lf %u\n", iter->id, phi, lambda, iter->offset);
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
    FILE *nodes = fopen(nodesFilename.c_str(), "r");
    FILE *edges = fopen(edgesFilename.c_str(), "r");
    FILE *output = fopen(outputFilename.c_str(), "w");

#ifdef DEBUG
    auto startNodes = std::chrono::steady_clock::now();
#endif
    std::vector<compactNode<uint>> V = proceedNodes(nodes);
    std::sort(V.begin(), V.end());
#ifdef DEBUG
    std::cout << "Nodes done in " << since(startNodes).count() / 1000 << " sec" << std::endl;
#endif


#ifdef DEBUG
    auto startEdges = std::chrono::steady_clock::now();
#endif
    std::vector<edge> E = proceedEdges(edges);
#ifdef DEBUG
    std::cout << "Edges read in " << since(startEdges).count() / 1000 << " sec" << std::endl;
#endif
    std::sort(E.begin(), E.end(), [](const edge &lhs, const edge &rhs) -> bool {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });
#ifdef DEBUG
    std::cout << "Edges done in " << since(startEdges).count() / 1000 << " sec" << std::endl;
#endif
    setOffset(output, V, E);
    for (int i = (int) V.size() - 2; i >= 0; --i) {
        V[i].offset = std::min(V[i].offset, V[i + 1].offset);
    }

    // Clear output file for correct offset
    fclose(output);
    output = fopen(outputFilename.c_str(), "w");

    fprintf(output, "%u\n", (uint) V.size());
    printNodes(output, nodes, V);
    printEdges(output, E, V);

    fclose(nodes);
    fclose(edges);
    fclose(output);
}

#endif //KP_DA_UTILS_H