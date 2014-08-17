#include "../src/arbpp.hpp"

#define BOOST_TEST_MODULE arb_test
#include <boost/test/unit_test.hpp>

// NOTE: this may change in the future to flint/flint.h, if Arb
// gets this change as well.
#include "flint.h"

using namespace arbpp;

BOOST_AUTO_TEST_CASE(arb_base_test)
{
    arb a0{20};
    a0 += 1;
    std::cout << a0 << '\n';
    arb a1{0.2};
    a1.set_precision(70);
    //a1.add_error(0.1f);
    std::cout << a1 << '\n';
    std::cout << (a0 + a1) << '\n';
    std::cout << (a0 + 6) << '\n';
    std::cout << (6 + a0) << '\n';
    std::cout << (6. + a0) << '\n';
    std::cout << (std::numeric_limits<long>::max() + a0) << '\n';
    a0.set_precision(100);
    std::cout << (std::numeric_limits<long>::max() + a0) << '\n';
    std::cout << (a0 + a1) << '\n';
    std::cout << cos(arb{0.0000001}) << '\n';
}

// Keep this for last, in order to have proper
// memory cleanup and make valgrind happy.
BOOST_AUTO_TEST_CASE(arb_cleanup)
{
    ::flint_cleanup();
}
