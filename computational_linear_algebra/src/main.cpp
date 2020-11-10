#include "matrix.hpp"
#include <iostream>
#include <string>
#include "les_solve.hpp"

using namespace lin_algebra;
using std::cin; using std::cout;


int main(){
    matrix<double> A(0, 3, 3), B(0, 3, 1), C(0, 3, 1);
    freopen("input.txt", "r", stdin);
    cin >> A;
    cin >> B;
    cin >> C;
    
    cout << iteration_solve<double>(A, B, 0.001);
}