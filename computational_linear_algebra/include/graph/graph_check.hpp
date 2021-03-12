#ifndef ALGEBRA_GRAPH_CHECK_HPP
#define ALGEBRA_GRAPH_CHECK_HPP

#include "../matrix/matrix.hpp"
#include "../eigenvalues/eigenvalues.hpp"
#include "../eigenvalues/tridiagonalization.hpp"
namespace lin_algebra {
    template<typename T>
    class graph : public real_and_complex_restriction<T> {
    public:
        graph(matrix<T> const &A_, long double e) : gr_(A_), eps(e), V(A_.get_rows()) {};
        ~graph() = default;

        std::vector<T> get_spectre() {
            if (!spectre_calced) {
                calc_spec();
            }

            return spectre_;
        }

        std::vector<int> get_degs() {
            if (!degs_calced) {
                calc_degs();
            }

            return degs_;
        }

        bool operator==(graph<T> &g) {
            if (V != g.V) return false;
            if (!degs_calced) calc_degs();
            if (!spectre_calced) calc_spec();

            std::vector<T> g_spec = g.get_spectre();
            std::vector<int> g_degs = g.get_degs();

            if (g_spec.size() != spectre_.size() || g_degs.size() != degs_.size()) return false;

            for(int i = 0; i < degs_.size(); i++) {
                if (degs_[i] != g_degs[i]) return false;
            }

            for(int i = 0; i < spectre_.size(); i++) {
                if (std::abs(g_spec[i] - spectre_[i]) > eps * 10) return false;
            }

            return true;
        }

        bool operator!=(graph<T> &g) {
            return !(*this == g);
        }
    private:
        void calc_spec() {
            triagonalize<T> tr(gr_);
            eigenvalues_solve<T> eg(tr.start().first);
            spectre_ = eg.QR_solve_shift(eps).first;
            std::sort(spectre_.begin(), spectre_.end());
            std::reverse(spectre_.begin(), spectre_.end());
            spectre_calced = true;
        }

        void calc_degs() {
            matrix<T> tmp(gr_ * gr_);
            for(int i = 0; i < tmp.get_rows(); i++) {
                degs_.push_back(tmp[i][i]);
            }
            std::sort(degs_.begin(), degs_.end());
            degs_calced = true;
        }

        bool spectre_calced = false, degs_calced = false;
        std::vector<T> spectre_;
        std::vector<int> degs_;
        matrix<T> gr_;
        long double eps;
        int V;
    };

    class graph_gen {
    public:
        graph_gen(int n_, long double e);

        long double ret_alpha();
    private:
        matrix<long double> gr;
        int d, n;
        long double eps;
    };

    class graph_gen_prime {
    public:
        graph_gen_prime(int p_, long double e);

        long double ret_alpha();
    private:
        matrix<long double> gr;
        int d, p;
        long double eps;
    };
}

#endif //ALGEBRA_GRAPH_CHECK_HPP
