#include "matrix.hpp"
#include <iostream>
#include <string>
#include "les_solve.hpp"
#include "qr_algo.hpp"

using namespace lin_algebra;
using std::cin; using std::cout;


int main(){
    matrix<long double> A(0, 3, 3), B(0, 3, 3), C(0, 3, 1);
    givens_matrix<long double> G(1, 2, 0.707, 0.707, 3);
    freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);
    cin >> A;
    cin >> B;
    cin >> C;
    QR_algo<long double> qr(A);
    auto Q = qr.givens_algo();
    LES_solver<long double> b(A, C);
    std::cout << b.gauss_seidel(0.000000001);
}