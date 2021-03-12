#include "../include/tests/test_class.hpp"
#include <iostream>
int main() {
    algebra_test a;
    a.run_all_tests();
    a.show_final_result();
}
