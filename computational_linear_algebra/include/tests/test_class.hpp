#ifndef ALGEBRA_TEST_CLASS_HPP
#define ALGEBRA_TEST_CLASS_HPP
#include <vector>
#include <string>

class test {
protected:
    static int failed_num; // количество тестов, которые сломались
    static int total_num;  // общее количество тестов

public:
    virtual void show_final_result() = 0;

    virtual void run_all_tests() = 0;
};


class algebra_test : public test {
private:
    std::vector<std::string> failed_list;

    void matrix_implementation_test_1();
    void matrix_implementation_test_2();
    void matrix_implementation_test_3();
    void matrix_implementation_test_4();

    void les_implementation_test_1();
    void les_implementation_test_2();
    void les_implementation_test_3();

    void qr_implementation_test_1();
    void qr_implementation_test_2();
    void qr_implementation_test_3();

    void eigenvalues_search_test_1();
    void eigenvalues_search_test_2();
    void eigenvalues_search_test_3();
    void eigenvalues_search_test_4();

    void graph_isom_test_1();
    void graph_isom_test_2();
    void graph_isom_test_3();
public:
    void run_all_tests();
    void show_final_result();
};


#endif //ALGEBRA_TEST_CLASS_HPP
