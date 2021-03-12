#include "../include/tests/test_class.hpp"
#include <iostream>
#include <fstream>
#include "../include/lin_algebra.hpp"
#include "../include/graph/graph_check.hpp"
using namespace lin_algebra;
using std::cout;
int test::total_num = 0;
int test::failed_num = 0;
long double epsilon = 1e-10;

void algebra_test::matrix_implementation_test_1() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/matrix_test1.txt");
    matrix<double> A1(0, 3, 3), A2(0, 3, 3), A3(0, 3, 3);
    in >> A1 >> A2 >> A3;

    if (A3 != A1 + A2) {
        failed_list.emplace_back("matrix_implementation_test_1");
        failed_num++;
    }
    total_num++;
}

void algebra_test::matrix_implementation_test_2() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/matrix_test2.txt");
    matrix<double> A1(0, 3, 3), A2(0, 3, 1), A3(0, 3, 1);
    in >> A1 >> A2 >> A3;

    if (A3 != A1 * A2) {
        failed_list.emplace_back("matrix_implementation_test_2");
        failed_num++;
    }
    total_num++;
}

void algebra_test::matrix_implementation_test_3() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/matrix_test3.txt");
    matrix<double> A1(0, 3, 4), A2(0, 4, 4), A3(0, 3, 4);
    in >> A1 >> A2 >> A3;

    if (A3 != A1 * A2) {
        failed_list.emplace_back("matrix_implementation_test_3");
        failed_num++;
    }
    total_num++;
}

void algebra_test::matrix_implementation_test_4() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/matrix_test4.txt");
    matrix<long double> A1(0, 4, 4);
    in >> A1;

    if (std::abs(A1.normilize().norm() - 1) > epsilon) {
        failed_list.emplace_back("matrix_implementation_test_4");
        failed_num++;
    }
    total_num++;
}

void algebra_test::les_implementation_test_1() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/les1.txt");
    matrix<double> A1(0, 2, 2), x(0, 2, 1), B(0, 2, 1);
    in >> A1;
    in >> x;
    in >> B;
    LES_solver<double> ls(A1, B);
    //cout << A1 << x << B;
    if (x != ls.gauss_seidel(epsilon, matrix<double>(0,1,1), true)) {
        failed_list.emplace_back("les_implementation_test_1");
        failed_num++;
    }
    total_num++;
}

void algebra_test::les_implementation_test_2() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/les2.txt");
    matrix<double> A1(0, 3, 3), x(0, 3, 1), B(0, 3, 1);
    in >> A1;
    in >> x;
    in >> B;
    LES_solver<double> ls(A1, B);

    if (x != ls.iteration_solve(epsilon, matrix<double>(0,1,1), true)) {
        failed_list.emplace_back("les_implementation_test_2");
        failed_num++;
    }
    total_num++;
}

void algebra_test::les_implementation_test_3() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/les3.txt");
    matrix<long double> A1(0, 4, 4), x(0, 4, 1), B(0, 4, 1);
    in >> A1;
    in >> x;
    in >> B;
    LES_solver<long double> ls(A1, B);
    if (x != ls.iteration_solve(1e-16, matrix<long double>(0,1,1), true)) {
        failed_list.emplace_back("les_implementation_test_3");
        failed_num++;
    }
    total_num++;
}

void algebra_test::qr_implementation_test_1() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/qr12.txt");
    matrix<double> A1(0, 4, 4);
    in >> A1;
    QR_algo<double> qr(A1);
    auto ans = qr.givens_algo();
    if (A1 != ans.first * ans.second) {
        failed_list.emplace_back("qr_implementation_test_1");
        failed_num++;
    }
    total_num++;
}

void algebra_test::qr_implementation_test_2() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/qr12.txt");
    matrix<double> A1(0, 4, 4);
    in >> A1;
    QR_algo<double> qr(A1);
    auto ans = qr.housholder_algo();

    if (A1 != ans.first * ans.second) {
        failed_list.emplace_back("qr_implementation_test_2");
        failed_num++;
    }
    total_num++;
}

void algebra_test::qr_implementation_test_3() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/qr3.txt");
    matrix<double> A1(0, 5, 5), Q(5);
    in >> A1;
    QR_algo<double> qr(A1);
    auto ans = qr.tridiag_algo();

    std::reverse(ans.first.begin(), ans.first.end());
    for(auto &e: ans.first) {
        Q = e.transpose() * Q;
    }

    if (A1 != Q * ans.second) {
        failed_list.emplace_back("qr_implementation_test_3");
        failed_num++;
    }
    total_num++;
}

void algebra_test::eigenvalues_search_test_1() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/eig1.txt");
    matrix<std::complex<long double>> A1(0, 4, 4);
    long double lambda;
    in >> A1 >> lambda;
    eigenvalues_solve<std::complex<long double>> eig(A1);

    if (std::abs(lambda - eig.iteration_solve(1e-8, matrix<std::complex<long double>>(1), true).first.real()) > 1e-4) {
        failed_list.emplace_back("eigenvalues_search_test_1");
        failed_num++;
    }
    total_num++;
}

void algebra_test::eigenvalues_search_test_2() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/eig2.txt");
    matrix<double> A1(0, 100, 100);
    double lambda1, lambda2;
    in >> A1 >> lambda1 >> lambda2;

    eigenvalues_solve<double> eig(A1);
    auto tmp = eig.QR_solve(1e-5);

    if (std::abs(lambda1 - tmp.first[0]) > 1e-4 || std::abs(lambda2 - tmp.first[1]) > 1e-4) {
        failed_list.emplace_back("eigenvalues_search_test_2");
        failed_num++;
    }
    total_num++;
}

void algebra_test::eigenvalues_search_test_3() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/eig3.txt");
    matrix<double> A1(0, 100, 100);
    double lambda1, lambda2;
    in >> A1 >> lambda1 >> lambda2;

    triagonalize<double> tr(A1);
    eigenvalues_solve<double> eig(tr.start().first);
    auto tmp = eig.QR_solve(1e-5);

    if (std::abs(lambda1 - tmp.first[0]) > 1e-4 || std::abs(lambda2 - tmp.first[1]) > 1e-4) {
        failed_list.emplace_back("eigenvalues_search_test_3");
        failed_num++;
    }
    total_num++;
}

void algebra_test::eigenvalues_search_test_4() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/eig4.txt");
    matrix<long double> A(0,10,10), D(0, 10, 10);
    std::vector<long double> lambdas(10);
    in >> A;
    for(int i = 0; i < 10; i++) {
        in >> lambdas[i];
    }

    bool check_lamd = true;
    triagonalize<long double> tr(A);
    matrix<long double> Q = tr.start().second;
    A = tr.start().first;
    eigenvalues_solve<long double> eig(A);
    auto tmp = eig.QR_solve_shift(1e-10);

    for (int i = 0; i < 10; i++) {
        if (abs(lambdas[i] - tmp.first[i]) > 1e-4) {
            check_lamd = false;
        }
        D[i][i] = tmp.first[i];
    }

    if (tmp.second * D * tmp.second.transpose() != A || check_lamd) {
        failed_list.emplace_back("eigenvalues_search_test_4");
        failed_num++;
    }
    total_num++;
}

void algebra_test::graph_isom_test_1() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/graph1.txt");
    matrix<long double> A1(0, 4, 4), A2(0,4,4);
    in >> A1 >> A2;
    graph<long double> gr(A1, 1e-5), gr2(A2, 1e-5);


    if (gr != gr2) {
        failed_list.emplace_back("graph_isom_test_1");
        failed_num++;
    }
    total_num++;
}

void algebra_test::graph_isom_test_2() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/graph2.txt");
    matrix<long double> A1(0, 4, 4), A2(0,4,4);
    in >> A1 >> A2;
    graph<long double> gr(A1, 1e-5), gr2(A2, 1e-5);

    if (gr == gr2) {
        failed_list.emplace_back("graph_isom_test_2");
        failed_num++;
    }
    total_num++;
}

void algebra_test::graph_isom_test_3() {
    std::ifstream in("/home/bagrorg/CLionProjects/algebra/tests/bank/graph3.txt");
    matrix<long double> A1(0, 5, 5), A2(0,5,5);
    in >> A1 >> A2;
    graph<long double> gr(A1, 1e-5), gr2(A2, 1e-5);

    if (gr != gr2) {
        failed_list.emplace_back("graph_isom_test_3");
        failed_num++;
    }
    total_num++;
}

void algebra_test::run_all_tests() {
    matrix_implementation_test_1();
    matrix_implementation_test_2();
    matrix_implementation_test_3();
    matrix_implementation_test_4();

    les_implementation_test_1();
    les_implementation_test_2();
    les_implementation_test_3();

    qr_implementation_test_1();
    qr_implementation_test_2();
    qr_implementation_test_3();

    eigenvalues_search_test_1();
    eigenvalues_search_test_2();
    eigenvalues_search_test_3();
    eigenvalues_search_test_4();

    graph_isom_test_1();
    graph_isom_test_2();
    graph_isom_test_3();
}

void algebra_test::show_final_result() {
    std::cout << "Total number of tests : " << total_num << "\n";
    std::cout << "Number of failed tests : " << failed_num << "\n";
    std::cout << "Failed tests : " << failed_num << "\n";
    for(auto &e: failed_list) {
        std::cout << "-- " << e << "\n";
    }
}

















