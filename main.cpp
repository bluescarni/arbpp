#include "arbpp.hpp"

using namespace arbpp;

int main()
{
    arb a0{20};
    std::cout << a0 << '\n';
    arb a1{0.2};
    a1.set_precision(70);
    //a1.add_error(0.1f);
    std::cout << a1 << '\n';
    std::cout << (a0 + a1) << '\n';
    a0.set_precision(100);
    std::cout << (a0 + a1) << '\n';
    std::cout << cos(arb{0.0000001}) << '\n';
    ::flint_cleanup();
}
