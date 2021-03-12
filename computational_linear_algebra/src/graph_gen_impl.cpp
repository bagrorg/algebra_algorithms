#include "../include/graph/graph_check.hpp"

namespace lin_algebra {
    long binpow (long a, long n, long m) {
        long res = 1;
        while (n) {
            if (n & 1) {
                res *= a;
                res %= m;
            }
            a *= (a % m);
            a %= m;
            n >>= 1;
        }
        return res % m;
    }

    graph_gen::graph_gen(int n_, long double e) {
        n = n_;
        gr = matrix<long double>(0, n_ * n_, n_ * n_);
        for (int x = 0; x < n; x++) {
            for(int y = 0; y < n; y++) {
                std::vector<std::pair<int, int>> neighbours = {
                        {x + 2 * y, y},
                        {x - 2 * y, y},
                        {x + 2 * y + 1, y},
                        {x - 2 * y - 1, y},
                        {x, y + 2 * x},
                        {x, y - 2 * x},
                        {x, y + 2 * x + 1},
                        {x, y - 2 * x - 1}
                };
                for (auto &e: neighbours) {
                    int nx = e.first, ny = e.second;
                    int v = ((x + n * n) % n) * n + (y + n * n) % n, u = ((nx + n * n) % n) * n + (ny + n * n) % n;
                    gr[v][u] += 1;
                }
            }
        }

        d = 0;
        for (int i = 0; i < n * n; i++) {
            d += gr[i][0];
        }
        eps = e;
    }

    long double graph_gen::ret_alpha() {
        graph<long double> tmp(gr, eps);
        auto spec = tmp.get_spectre();
        if (spec.size() == 0) return -1;
        long double alpha = std::max(std::abs(spec[1]), std::abs(spec[spec.size() - 1])) / d;
        return alpha;
    }

    graph_gen_prime::graph_gen_prime(int p_, long double e) {
        p = p_;
        gr = matrix<long double>(0, p_ + 1, p_ + 1);
        for (int x = 1; x < p; x++) {
            gr[x][(x + 1) % p] += 1;
            gr[x][(x - 1) % p] += 1;
            gr[x][binpow(x, p-2, p)] += 1;
        }
        gr[p][0] += 1;
        gr[p][p] += 1;
        gr[p][p] += 1;

        gr[0][p] += 1;
        gr[0][p-1] += 1;
        gr[0][1] += 1;

        d = 0;
        for (int i = 0; i < p + 1; i++) {
            d += gr[i][0];
        }
        eps = e;
    }

    long double graph_gen_prime::ret_alpha() {
        graph<long double> tmp(gr, eps);
        auto spec = tmp.get_spectre();
        if (spec.size() == 0) return -1;
        long double alpha = std::max(std::abs(spec[1]), std::abs(spec[spec.size() - 1])) / d;
        return alpha;
    }
}
