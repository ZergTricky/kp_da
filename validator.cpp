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
const ld eps = 1e-6;
const ld PI = acos(-1);

const double R = 6371;

inline double rad(double angle) {
    return PI / 180.0 * angle;
}

double dist(const pair<double, double> &a, const pair<double, double> &b) {
    double cosd = sin(a.first) * sin(b.first) + cos(a.first) * cos(b.first) * cos(a.second - b.second);
    assert(fabs(cosd) < 1 + eps);
    double d = acos(cosd);
    return R * d;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    if (argc < 5) {
        cout << "Not enough params" << "\n";
        return 0;
    }
    ifstream nodes(argv[1]);
    ifstream edges(argv[2]);
    ifstream queries(argv[3]);
    ifstream res(argv[4]);
    map<ll, pair<double, double>> data;
    ll id;
    while (nodes >> id) {
        data[id];
        nodes >> data[id].first >> data[id].second;
        data[id].first = rad(data[id].first);
        data[id].second = rad(data[id].second);
    }
    map<ll, vector<ll>> g;
    ll k;
    while (edges >> k) {
        vector<ll> path;
        path.resize(k);
        for (auto &x: path) {
            edges >> x;
        }
        for (ll i = 1; i < k; ++i) {
            g[path[i]].push_back(path[i - 1]);
            g[path[i - 1]].push_back(path[i]);
        }
    }
    vector<pl> q;
    ll x, y;
    while (queries >> x >> y) {
        q.push_back({x, y});
    }
    ll cur = 0;
    ld d;
    while (res >> d) {
        ll n;
        res >> n;
        vector<ll> way(n);
        for (auto &x: way) res >> x;

        ld real = 0.0;
        for (ll i = 1; i < n; ++i) {
            real += dist(data[way[i - 1]], data[way[i]]);
            bool ok = false;
            for (auto &to: g[way[i]]) {
                if (to == way[i - 1])ok = true;
            }
            if (!ok) {
                cout << "BAD: " << q[cur].first << " " << q[cur].second << "\n";
                exit(0);
            }
        }
        if (fabs(d + 1) > eps) {
            assert(fabs(real - d) < eps);
        }

        map<ll, double> D;
        set<pair<double, ll>> s;
        D[q[cur].first] = 0;
        s.insert({D[q[cur].first], q[cur].first});
        while (!s.empty()) {
            auto it = *s.begin();
            s.erase(s.begin());
            ll v = it.second;
            for (auto to: g[v]) {
                double w = dist(data[v], data[to]);
                if (!D.contains(to) || D[to] > D[v] + w + eps) {
                    s.erase({D[to], to});
                    D[to] = D[v] + w;
                    s.insert({D[to], to});
                }
            }
        }
        if (!D.contains(q[cur].second)) {
            assert(fabs(d + 1) < eps);
        } else {
            std::cout << d << " " << D[q[cur].second] << std::endl;
            assert(fabs(d - D[q[cur].second]) < eps);
        }

        ++cur;
    }
    cout << "VALIDATION OK" << endl;
    return 0;
}