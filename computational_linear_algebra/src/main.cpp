#include "matrix.hpp"
#include <iostream>
#include <string>
#include "les_solve.hpp"

using namespace lin_algebra;
using std::cin; using std::cout;


int main(){
    matrix<float> A(0, 3, 3), B(0, 3, 1);
    freopen("input.txt", "r", stdin);
    cin >> A;
    cin >> B;
    matrix<float> C(A);
    LES_solver<float> a(A, B);
    cout << a.gauss_seidel(0.000001);
}