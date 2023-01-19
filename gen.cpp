#include <bits/stdc++.h>

using namespace std;

#define flush cout.flush

using ll = long long;
using ull = unsigned long long;
using ld = long double;
using pl = pair<ll, ll>;
const ll INF = 1e9 + 7;
const ll mod = 1e9 + 7;
const ll mod2 = 998244353;
const ld eps = 1e-9;
const ld PI = acos(-1);

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    if (argc < 4) {
        std::cout << "Not enough params" << std::endl;
        return 0;
    }
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    int n = stoi(argv[1]);
    int m = stoi(argv[2]);
    int q = stoi(argv[3]);

    ofstream nodes("n_gen.txt");
    ofstream edges("e_gen.txt");
    ofstream queries("q_gen.txt");

    vector<ll> p(n);
    iota(p.begin(), p.end(), 0);
    shuffle(p.begin(), p.end(), rng);

    for (int i = 0; i < n; ++i) {
        nodes << p[i] + 1 << " " << rng() % 90 << " " << rng() % 90 << "\n";
    }
    for (int i = 0; i < m; ++i) {
        int k = rng() % 4 + 2;
        edges << k << " ";
        for (int j = 0; j < k; ++j) {
            edges << rng() % n + 1 << " ";
        }
        edges << "\n";
    }
    for (int i = 0; i < q; ++i) {
        queries << rng() % n + 1 << " " << rng() % n + 1 << "\n";
    }
    return 0;
}