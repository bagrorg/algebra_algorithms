#ifndef EIGENVALUES_HPP
#define EIGENVALUES_HPP

#include "../matrix/matrix.hpp"
#include "../qr/qr_algo.hpp"

namespace lin_algebra {
    using std::abs, std::sqrt;
    //Если вещественные, то QR алгоритм
    template<typename T>
    class eigenvalues_solve : public real_and_complex_restriction<T> {
    public:
        //Матрица A
        eigenvalues_solve(const matrix <T> &A) : n_(A.get_rows()), A_(A) {};

        ~eigenvalues_solve() = default;

        //QR алгоритм
        //Передается точность
        std::pair<std::vector<T>, matrix <T>> QR_solve(long double eps) {
            matrix <T> Ai(A_), Qi(n_);
            int steps = 0;
            //Пока радиусы кругов не станут хорошими
            while (!Ai.gersh_small_circle(eps)) {
                //Если долго работаем, то выводим 0
                if (steps > STEP_LIMITS) return {{}, matrix<T>(0,1,1)};
                steps++;


                QR_algo <T> qr(Ai);
                //Делаем QR разложение
                auto tmp = qr.housholder_algo();
                //Меняем местами Q и R
                Ai = tmp.second * tmp.first;
                //Домножаем Qi на Q
                Qi = Qi * tmp.first;
            }
            //Вектор собственных чисел
            std::vector<T> ret;
            for(int i = 0; i < Ai.get_rows(); i++) {
                ret.push_back(Ai[i][i]);
            }

            return {ret, Qi};
        }
        //QR алгоритм для трехдиагональных матриц
        //Передается точность
        std::pair<std::vector<T>, matrix <T>> QR_solve_triag(long double eps) {
            matrix <T> Ai(A_), Qi(n_);
            int steps = 0;
            //Пока радиусы кругов не станут хорошими
            while (!Ai.gersh_small_circle(eps)) {
                //Если долго работаем, то выводим 0
                if (steps > STEP_LIMITS) return {{}, matrix<T>(0,1,1)};
                steps++;
                QR_algo <T> qr(Ai);
                //Делаем QR разложение
                std::pair<std::vector<givens_matrix<T>>, matrix<T>> tmp = qr.tridiag_algo();
                //Домножаем n раз за n (так как гивенсом умеем быстро) - n^2 действий
                for(givens_matrix<T> &G : tmp.first) {
                    tmp.second = G.mult_from_right(tmp.second);
                    Qi = G * Qi;
                }
                //Обновляем матрицу
                Ai = tmp.second;
            }
            //Список сч
            std::vector<T> ret;
            for(int i = 0; i < Ai.get_rows(); i++) {
                ret.push_back(Ai[i][i]);
            }

            return {ret, Qi.transpose()};
        }
        //QR алгоритм с использованием сдвига
        //Передается точность
        std::pair<std::vector<T>, matrix<T>> QR_solve_shift(long double eps) {
            int steps = 0;
            matrix <T> Ai(A_), Qi(n_);
            std::vector<T> ret;
            //Пока не дойдем до подматрицы размера 1х1
            while (Ai.get_rows() > 1) {
                //Если долго работаем, то выводим 0
                if (steps > STEP_LIMITS) return {{}, matrix<T>(0,1,1)};
                steps++;


                //Список из элементов правого нижнего квадрата
                std::vector<T> cf;
                for(int i = Ai.get_rows() - 2; i < Ai.get_rows(); i++) {
                    for(int j = Ai.get_cols() - 2; j < Ai.get_cols(); j++) {
                        cf.push_back(Ai[i][j]);
                    }
                }
                //Ищем сч правого нижнего квадрата через квадратное уравнение
                T D = (cf[0] - cf[3]) * (cf[0] - cf[3]) + 4 * cf[1] * cf[2];
                T x1 = (cf[0] + cf[3] + sqrt(D)) / 2, x2 = (cf[0] + cf[3] - sqrt(D)) / 2;

                //Приближение -- берем ближайшее к a_n,n
                T vk;
                if (D < 0) vk = Ai[Ai.get_rows() - 1][Ai.get_rows() - 1];
                else if (abs(x1 - cf[3]) < abs(x2 - cf[3])) {
                    vk = x1;
                } else {
                    vk = x2;
                }
                //Единичная матрица
                matrix<T> En(Ai.get_rows());

                //Заводим куар разложение
                QR_algo <T> qr(Ai - (En * vk));
                //Делаем QR разложение
                std::pair<std::vector<givens_matrix<T>>, matrix<T>> tmp = qr.tridiag_algo();

                for(givens_matrix<T> &G : tmp.first) {
                    tmp.second = G.mult_from_right(tmp.second); //Ai = Ai * G
                    Qi = G * Qi;
                }
                Ai = tmp.second;
                //Сдвигаем обратно
                Ai += (En * vk);
                //Проверяем, что столбец и строка последние достаточно маленькие
                T bot = 0, right = 0;
                for(int i = 0; i < Ai.get_rows() - 1; i++) {
                    bot += abs(Ai[Ai.get_rows() - 1][i]);
                    right += abs(Ai[i][Ai.get_cols() - 1]);
                }
                //Если достаточно маленькие, то переходим к меньшей подматрице
                if (bot < eps && right < eps) {
                    ret.push_back(Ai[Ai.get_rows() - 1][Ai.get_cols() - 1]);
                    Ai = Ai.delete_last_row_col();
                }
            }
            ret.push_back(Ai[0][0]);
            std::reverse(ret.begin(), ret.end());
            return {ret, Qi.transpose()};
        }
    private:
        std::size_t n_;
        matrix <T> A_;
        int bad_count = 0;
        const int STEP_LIMITS = 100000;
    };

    using comp = std::complex<long double>;

    //Если комплексный шаблон, то используем метод итераций
    template<>
    class eigenvalues_solve<comp> : public real_and_complex_restriction<comp> {
    public:
        eigenvalues_solve(const matrix <comp> &A) : n_(A.get_rows()), A_(A) {};

        ~eigenvalues_solve() = default;

        //Метод итераций (точность, приближение, надо ли генерировать случайное приближение)
        std::pair<comp, matrix < comp>> iteration_solve(long double eps, matrix<comp> x0, bool generate_random) {
            //Генерируем случайный комплексный вектор
            if (generate_random) x0 = matrix < comp > ::random_generate(0.5, n_, 1);
            //СЧ которое ищем
            comp lambda = (x0.transpose() * A_ * x0)[0][0];
            //Пока норма разницы не будет достаточно маленькой
            while ((A_ * x0 - (x0 * lambda)).norm() >= eps) {
                //Берем новый вектор нормализованным//Если долго работаем, то выводим 0
                matrix < comp > x_new = (A_ * x0).normilize();
                //Если расходимся уже долго, то выходим из итераций
                if (abs((A_ * x0 - x0 * lambda).norm()) <=
                    abs((A_ * x_new - x_new * lambda).norm()) || abs((A_ * x0 - x0 * lambda).norm()) -
                                                                      abs((A_ * x_new - x_new * lambda).norm()) <=
                                                                      1e-7) bad_count++;
                else bad_count = 0;
                //Если долго расходимся, то возвращаем 0
                if (bad_count == FAILURE_LIMIT) return {0, matrix < comp > (0, 1, 1)};

                x0 = x_new;
                lambda = (x0.transpose() * A_ * x0)[0][0];
            }

            return {lambda, x0};
        }

    private:
        std::size_t n_;
        matrix <comp> A_;
        int bad_count = 0;
        const int FAILURE_LIMIT = 20;
    };
}


#endif //ALGEBRA_EIGENVALUES_HPP
