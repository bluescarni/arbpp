#include "../src/arbpp.hpp"

#define BOOST_TEST_MODULE arb_test
#include <boost/test/unit_test.hpp>

#include "flint/flint.h"

using namespace arbpp;

BOOST_AUTO_TEST_CASE(arb_base_test)
{
    arb a;
}

// Keep this for last, in order to have proper
// memory cleanup and make valgrind happy.
BOOST_AUTO_TEST_CASE(arb_cleanup)
{
    ::flint_cleanup();
}
