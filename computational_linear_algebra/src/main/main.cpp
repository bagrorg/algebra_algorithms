#include <string>
#include "../../include/lin_algebra.hpp"
#include "../../include/graph/graph_check.hpp"
using namespace lin_algebra;
using std::cin; using std::cout;
using std::string;

int main(){
    cout << "Введите txt файл ввода >> ";
    string file;
    cin >> file;
    string command;
    while (true) {
        cout << std::fixed << "Введите номер Задания или exit\n>> ";
        cin >> command;
        if (command == "exit") {
            break;
        } else if (command == "1") {
            cout << "Введите квадратную матрицу A (такую что x = Ax + b), вектор b, eps - точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n), b(0, n, 1);
                long double eps;
                std::ifstream in(file);

                in >> A >> b >> eps;

                LES_solver<long double> les(A, b);
                matrix<long double> ans = les.iteration_solve(eps, matrix<long double>(1), true);
                if (ans == matrix<long double>(0,1,1)) {
                    cout << "Метод разошелся.. :(\n\n";
                } else {
                    cout << ans << '\n';
                }
                in.close();
            }
        } else if (command == "2") {
            cout << "Введите квадратную матрицу A (такую что Ax = b), вектор b, eps - точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n), b(0, n, 1);
                long double eps;
                std::ifstream in(file);

                in >> A >> b >> eps;
                LES_solver<long double> les(A, b);
                matrix<long double> ans = les.gauss_seidel(eps, matrix<long double>(1), true);
                if (ans == matrix<long double>(0,1,1)) {
                    cout << "Метод разошелся.. :(\n\n";
                } else {
                    cout << ans << '\n';
                }
                in.close();
            }
        } else if (command == "3") {
            cout
                    << "Введите квадратную матрицу A, i j - коэффициенты, c и s такие, что c^2 + s^2 = 1 в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n, i, j;
                cin >> n;
                matrix<long double> A(0, n, n), b(0, n, 1);
                long double c, s;

                std::ifstream in(file);
                in >> A >> i >> j >> c >> s;

                givens_matrix<long double> G(i, j, s, c, n);
                cout << G * A << '\n';
                in.close();
            }
        } else if (command == "4") {
            cout << "Введите квадратную матрицу A в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                std::ifstream in(file);

                in >> A;
                QR_algo<long double> qr(A);
                auto ans = qr.givens_algo();
                cout << "Q = \n" << ans.first << "\nR = \n" << ans.second << "\n";
                in.close();
            }
        } else if (command == "5") {
            cout << "Введите квадратную матрицу A и вектор v в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n), v(0, n, 1);
                std::ifstream in(file);

                in >> A >> v;
                householder_matrix<long double> H(n, v);
                auto ans = H * A;
                cout << ans << "\n";
                in.close();
            }
        } else if (command == "6") {
            cout << "Введите квадратную матрицу A в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                std::ifstream in(file);

                in >> A;
                QR_algo<long double> qr(A);
                auto ans = qr.housholder_algo();
                cout << "Q = \n" << ans.first << "\nR = \n" << ans.second << "\n";
                in.close();
            }
        } else if (command == "7") {
            cout << "Введите квадратную симметричную матрицу A, вектор x0 и точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<std::complex<long double>> A(0, n, n), x0(0, n, 1);
                long double eps;
                std::ifstream in(file);

                in >> A >> x0 >> eps;
                eigenvalues_solve<std::complex<long double>> eig(A);
                auto ans = eig.iteration_solve(eps, x0, false);
                if (ans.first == (comp) 0) {
                    cout << "Метод разошелся.. :(\n\n";
                    break;
                }
                cout << "lambda = " << ans.first << "\n" << "v = \n" << ans.second << "\n";
                in.close();
            }
        } else if (command == "8") {
            cout << "Введите квадратную симметричную матрицу A и точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                long double eps;
                std::ifstream in(file);

                in >> A >> eps;
                eigenvalues_solve<long double> eig(A);
                auto ans = eig.QR_solve(eps);
                if (ans.first.empty()) {
                    cout << "Метод разошелся.. :(\n\n";
                }
                cout << "lambdas : ";
                for (auto e: ans.first) {
                    cout << e << ' ';
                }
                cout << '\n';
                cout << "Q^(k) = \n" << ans.second << "\n\n";
                in.close();
            }
        } else if (command == "9") {
            cout << "Введите квадратную симметричную матрицу A в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                std::ifstream in(file);
                in >> A;
                triagonalize<long double> tr(A);
                auto ans = tr.start();
                cout << "A' = \n" << ans.first << "\nQ = \n" << ans.second << "\n\n";
                in.close();
            }
        } else if (command == "10") {
            cout << "Введите квадратную трехдиагональную матрицу A и точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                long double eps;
                std::ifstream in(file);

                in >> A >> eps;

                eigenvalues_solve<long double> eig(A);
                auto ans = eig.QR_solve_triag(eps);
                if (ans.first.empty()) {
                    cout << "Метод разошелся.. :(\n\n";
                }
                cout << "lambdas : ";
                for (auto e: ans.first) {
                    cout << e << ' ';
                }
                cout << '\n';
                cout << "Q^(k) = \n" << ans.second << "\n\n";
                in.close();
            }
        } else if (command == "11") {
            cout << "Введите квадратную трехдиагональную матрицу A и точность в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размер матрицы >> ";
                int n;
                cin >> n;
                matrix<long double> A(0, n, n);
                long double eps;
                std::ifstream in(file);

                in >> A >> eps;
                eigenvalues_solve<long double> eig(A);
                auto ans = eig.QR_solve_shift(eps);
                if (ans.first.empty()) {
                    cout << "Метод разошелся.. :(\n\n";
                }
                cout << "lambdas : ";
                for (auto e: ans.first) {
                    cout << e << ' ';
                }
                cout << '\n';
                cout << "Q^(k) = \n" << ans.second << "\n\n";
                in.close();
            }
        } else if (command == "12") {
            cout << "Введите матрицы смежности графов в файл `input.txt`\nВвод готов? (напишите Y, когда будете готовы, N чтобы отменить) >> ";
            cin >> command;
            if (command == "Y" || command == "y") {
                cout << "Введите размеры матриц >> ";
                int n, m;
                cin >> n >> m;
                matrix<long double> A(0, n, n), B(0, m, m);
                long double eps = 1e-5;
                std::ifstream in(file);

                in >> A >> B;
                graph<long double> gr1(A, eps), gr2(B, eps);
                if (gr1 == gr2) {
                    cout << "Изоморфны \n";
                } else {
                    cout << "Неизоморфны \n";
                }
                in.close();
            }
        } else if (command == "13") {
            cout << "Введите подпункт >> ";
            int pod;
            cin >> pod;
            switch (pod) {
                case 1: {
                    cout << "Введите n >> ";
                    int n;
                    cin >> n;
                    graph_gen g(n, 1e-4);
                    cout << "alpha = " << g.ret_alpha() << "\n\n";
                    break;
                }
                case 2: {
                    cout << "Введите p >> ";
                    int p;
                    cin >> p;
                    graph_gen_prime gp(p, 1e-4);
                    cout << "alpha = " << gp.ret_alpha() << "\n\n";
                    break;
                }
                default: {
                    cout << "Unknown command\n";
                    break;
                }
            }
        } else {
            cout << "Unknown command!\n\n";
        }
    }
}
