/***************************************************************************
 *   Copyright (C) 2014 by Francesco Biscani                               *
 *   bluescarni@gmail.com                                                  *
 *                                                                         *
 *   This file is part of Arbpp.                                           *
 *                                                                         *
 *   Arbpp is free software: you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Arbpp is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.       *
 ***************************************************************************/

#include "../src/arbpp.hpp"

#define BOOST_TEST_MODULE arb_test
#include <boost/test/unit_test.hpp>

#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "arb.h"
#include "arf.h"
#include "mag.h"
#include "flint/flint.h"

using namespace arbpp;

BOOST_AUTO_TEST_CASE(arb_ctor_assignment_test)
{
    // Some type-traits checks.
    BOOST_CHECK((!std::is_constructible<arb,long double>::value));
    BOOST_CHECK((!std::is_constructible<arb,long long>::value));
    BOOST_CHECK((std::is_constructible<arb,long>::value));
    BOOST_CHECK((std::is_constructible<arb,char>::value));
    BOOST_CHECK((std::is_constructible<arb,unsigned char>::value));
    BOOST_CHECK((std::is_assignable<arb &,unsigned char>::value));
    BOOST_CHECK((std::is_assignable<arb &,double>::value));
    BOOST_CHECK((!std::is_assignable<arb &,long double>::value));
    // Default ctor.
    arb a0;
    BOOST_CHECK(::arf_is_zero(arb_midref(a0.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a0.get_arb_t())));
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision());
    // Copy ctor.
    arb a1{a0};
    BOOST_CHECK(::arf_is_zero(arb_midref(a1.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a1.get_arb_t())));
    BOOST_CHECK_EQUAL(a1.get_precision(),arb::get_default_precision());
    a1 = arb{1};
    a1.set_precision(100);
    arb a2{a1};
    BOOST_CHECK(::arf_is_one(arb_midref(a2.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a2.get_arb_t())));
    BOOST_CHECK_EQUAL(a2.get_precision(),100);
    // Move ctor.
    arb a3{std::move(a2)};
    BOOST_CHECK(::arf_is_one(arb_midref(a3.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a3.get_arb_t())));
    BOOST_CHECK_EQUAL(a3.get_precision(),100);
    BOOST_CHECK(::arf_is_zero(arb_midref(a2.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a2.get_arb_t())));
    BOOST_CHECK_EQUAL(a2.get_precision(),arb::get_default_precision());
    // Generic ctor.
    arb a4{42};
    BOOST_CHECK_EQUAL(a4.get_midpoint(),42.);
    BOOST_CHECK_EQUAL(a4.get_radius(),0.);
    BOOST_CHECK_EQUAL(a4.get_radius(),-0.);
    BOOST_CHECK_EQUAL(a4.get_precision(),arb::get_default_precision());
    BOOST_CHECK_EQUAL(arb{-42}.get_midpoint(),-42.);
    BOOST_CHECK_EQUAL(arb{-42}.get_radius(),0);
    BOOST_CHECK_EQUAL(arb{-42}.get_precision(),arb::get_default_precision());
    BOOST_CHECK_EQUAL(arb{12u}.get_midpoint(),12.);
    BOOST_CHECK_EQUAL(arb{12ul}.get_radius(),0);
    BOOST_CHECK_EQUAL(arb{12ul}.get_precision(),arb::get_default_precision());
    BOOST_CHECK_EQUAL(arb{1.3}.get_midpoint(),1.3);
    BOOST_CHECK_EQUAL(arb{1.3}.get_radius(),0.);
    BOOST_CHECK_EQUAL(arb{1.3}.get_precision(),arb::get_default_precision());
    // Generic constructor from precision.
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() + 1}.get_midpoint()),-42.);
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() + 1}.get_radius()),0);
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() + 1}.get_precision()),arb::get_default_precision() + 1);
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() - 1}.get_midpoint()),-42.);
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() - 1}.get_radius()),0);
    BOOST_CHECK_EQUAL((arb{-42,arb::get_default_precision() - 1}.get_precision()),arb::get_default_precision() - 1);
    // Copy assignment.
    arb a5;
    a5 = a4;
    BOOST_CHECK_EQUAL(a5.get_midpoint(),42.);
    BOOST_CHECK_EQUAL(a5.get_radius(),0.);
    BOOST_CHECK_EQUAL(a5.get_precision(),arb::get_default_precision());
    a4.set_precision(100);
    a4.add_error(1);
    a5 = a4;
    BOOST_CHECK_EQUAL(a5.get_midpoint(),42.);
    // NOTE: it seems that operations on mag_t, used to represent the radius,
    // will not be exact, and that the results might be a few ulp away (but
    // always producing strict upper or lower bounds). Just check that it is
    // not zero.
    BOOST_CHECK(a5.get_radius() != 0);
    BOOST_CHECK_EQUAL(a5.get_precision(),100);
    // Move assignment.
    a4.set_precision(101);
    a5 = std::move(a4);
    BOOST_CHECK_EQUAL(a5.get_midpoint(),42.);
    BOOST_CHECK(a5.get_radius() != 0);
    BOOST_CHECK_EQUAL(a5.get_precision(),101);
    BOOST_CHECK_EQUAL(a4.get_midpoint(),42.);
    BOOST_CHECK(a4.get_radius() != 0);
    BOOST_CHECK_EQUAL(a4.get_precision(),100);
    // Generic assignment.
    a1 = arb{0.5};
    a1.set_precision(100);
    a1.add_error(.1);
    BOOST_CHECK(!::mag_is_zero(arb_radref(a1.get_arb_t())));
    a1 = 1;
    BOOST_CHECK(::arf_is_one(arb_midref(a1.get_arb_t())));
    BOOST_CHECK(::mag_is_zero(arb_radref(a1.get_arb_t())));
    BOOST_CHECK_EQUAL(a1.get_precision(),53);
}

BOOST_AUTO_TEST_CASE(arb_add_error_test)
{
    arb a0;
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    a0.add_error(.1);
    BOOST_CHECK(a0.get_radius() >= .1);
    if (std::numeric_limits<double>::has_infinity) {
        a0.add_error(std::numeric_limits<double>::infinity());
        BOOST_CHECK(a0.get_radius() == std::numeric_limits<double>::infinity());
    }
    if (std::numeric_limits<double>::has_quiet_NaN) {
        BOOST_CHECK_THROW(a0.add_error(std::numeric_limits<double>::quiet_NaN()),std::invalid_argument);
    }
    BOOST_CHECK_THROW(a0.add_error(-1.),std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(arb_precision_test)
{
    arb a0{1};
    a0.set_precision(30);
    BOOST_CHECK_EQUAL(a0.get_precision(),30);
    BOOST_CHECK_THROW(a0.set_precision(0),std::invalid_argument);
    BOOST_CHECK_THROW(a0.set_precision(-1),std::invalid_argument);
    BOOST_CHECK_EQUAL(a0.get_precision(),30);
    BOOST_CHECK_EQUAL(a0.get_midpoint(),1.);
}

BOOST_AUTO_TEST_CASE(arb_get_arb_t_test)
{
    arb a0{1};
    BOOST_CHECK(a0.get_arb_t() != nullptr);
    BOOST_CHECK(static_cast<const arb &>(a0).get_arb_t() != nullptr);
}

BOOST_AUTO_TEST_CASE(arb_swap_test)
{
    arb a0{1}, a1{100};
    a0.set_precision(30);
    a0.add_error(.4);
    a0.swap(a0);
    BOOST_CHECK_EQUAL(a0.get_precision(),30);
    BOOST_CHECK_EQUAL(a0.get_midpoint(),1.);
    BOOST_CHECK(a0.get_radius() >= .4);
    a0.swap(a1);
    BOOST_CHECK_EQUAL(a1.get_precision(),30);
    BOOST_CHECK_EQUAL(a1.get_midpoint(),1.);
    BOOST_CHECK(a1.get_radius() >= .4);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision());
    BOOST_CHECK_EQUAL(a0.get_midpoint(),100.);
    BOOST_CHECK(a0.get_radius() == 0.);
}

BOOST_AUTO_TEST_CASE(arb_stream_test)
{
    arb a0{123.456};
    a0.add_error(.5);
    std::ostringstream oss;
    oss << a0;
    // Just check that some output is produced.
    BOOST_CHECK(oss.str().size() != 0);
}

BOOST_AUTO_TEST_CASE(arb_arithmetic_test)
{
    // In-place addition.
    arb a0{1}, a1{2};
    a0 += a1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),3.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision());
    // Try with different precisions.
    a1.set_precision(arb::get_default_precision() + 10);
    a0 += a1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),5.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // With plain int and unsigned.
    a0 += 1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),6.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    a0 += 1u;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),7.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // Float and double.
    a0 += 1.f;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),8.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    a0 += 2.;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),10.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // Binary add.
    arb a2{3}, a3{-4};
    BOOST_CHECK_EQUAL((a2 + a3).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((a2 + a3).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + a3).get_precision(),arb::get_default_precision());
    // Different precisions.
    a2.set_precision(arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 + a3).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((a2 + a3).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + a3).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a3 + a2).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((a3 + a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((a3 + a2).get_precision(),arb::get_default_precision() + 20);
    // With int and unsigned.
    BOOST_CHECK_EQUAL((a2 + 1).get_midpoint(),4.);
    BOOST_CHECK_EQUAL((a2 + 1).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + 1).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((1 + a2).get_midpoint(),4.);
    BOOST_CHECK_EQUAL((1 + a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((1 + a2).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 + 2u).get_midpoint(),5.);
    BOOST_CHECK_EQUAL((a2 + 2u).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + 2u).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((2u + a2).get_midpoint(),5.);
    BOOST_CHECK_EQUAL((2u + a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((2u + a2).get_precision(),arb::get_default_precision() + 20);
    // With floating-point.
    BOOST_CHECK_EQUAL((a2 + 1.f).get_midpoint(),4.);
    BOOST_CHECK_EQUAL((a2 + 1.f).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + 1.f).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((1.f + a2).get_midpoint(),4.);
    BOOST_CHECK_EQUAL((1.f + a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((1.f + a2).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 + 2.).get_midpoint(),5.);
    BOOST_CHECK_EQUAL((a2 + 2.).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 + 2.).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((2. + a2).get_midpoint(),5.);
    BOOST_CHECK_EQUAL((2. + a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((2. + a2).get_precision(),arb::get_default_precision() + 20);
    // In-place subtraction.
    a0.set_precision(arb::get_default_precision());
    a0 = 1;
    a1.set_precision(arb::get_default_precision());
    a1 = 2;
    a0 -= a1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-1.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision());
    // Try with different precisions.
    a1.set_precision(arb::get_default_precision() + 10);
    a0 -= a1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-3.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // With plain int and unsigned.
    a0 -= 1;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-4.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    a0 -= 1u;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-5.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // Float and double.
    a0 -= 1.f;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-6.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    a0 -= 2.;
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-8.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(a0.get_precision(),arb::get_default_precision() + 10);
    // Binary sub.
    a2 = arb{3};
    a3 = arb{4};
    BOOST_CHECK_EQUAL((a2 - a3).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((a2 - a3).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - a3).get_precision(),arb::get_default_precision());
    // Different precisions.
    a2.set_precision(arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 - a3).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((a2 - a3).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - a3).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a3 - a2).get_midpoint(),1.);
    BOOST_CHECK_EQUAL((a3 - a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((a3 - a2).get_precision(),arb::get_default_precision() + 20);
    // With int and unsigned.
    BOOST_CHECK_EQUAL((a2 - 1).get_midpoint(),2.);
    BOOST_CHECK_EQUAL((a2 - 1).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - 1).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((1 - a2).get_midpoint(),-2.);
    BOOST_CHECK_EQUAL((1 - a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((1 - a2).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 - 2u).get_midpoint(),1.);
    BOOST_CHECK_EQUAL((a2 - 2u).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - 2u).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((2u - a2).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((2u - a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((2u - a2).get_precision(),arb::get_default_precision() + 20);
    // With floating-point.
    BOOST_CHECK_EQUAL((a2 - 1.f).get_midpoint(),2.);
    BOOST_CHECK_EQUAL((a2 - 1.f).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - 1.f).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((1.f - a2).get_midpoint(),-2.);
    BOOST_CHECK_EQUAL((1.f - a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((1.f - a2).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((a2 - 2.).get_midpoint(),1.);
    BOOST_CHECK_EQUAL((a2 - 2.).get_radius(),0.);
    BOOST_CHECK_EQUAL((a2 - 2.).get_precision(),arb::get_default_precision() + 20);
    BOOST_CHECK_EQUAL((2. - a2).get_midpoint(),-1.);
    BOOST_CHECK_EQUAL((2. - a2).get_radius(),0.);
    BOOST_CHECK_EQUAL((2. - a2).get_precision(),arb::get_default_precision() + 20);
}

BOOST_AUTO_TEST_CASE(arb_negate_test)
{
    arb a0;
    a0.negate();
    BOOST_CHECK_EQUAL(a0.get_midpoint(),0.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    a0 = 42;
    a0.negate();
    BOOST_CHECK_EQUAL(a0.get_midpoint(),-42.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    a0.negate();
    BOOST_CHECK_EQUAL(a0.get_midpoint(),42.);
    BOOST_CHECK_EQUAL(a0.get_radius(),0.);
    BOOST_CHECK_EQUAL(-a0.get_midpoint(),-42.);
    BOOST_CHECK_EQUAL(-a0.get_radius(),0.);
    BOOST_CHECK_EQUAL((-(-a0)).get_midpoint(),42.);
    BOOST_CHECK_EQUAL((-(-a0)).get_radius(),0.);
}

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
    a0.add_error(1./0.);
    std::cout << a0 << '\n';
    a0.get_midpoint();
}

// Keep this for last, in order to have proper
// memory cleanup and make valgrind happy.
BOOST_AUTO_TEST_CASE(arb_cleanup)
{
    ::flint_cleanup();
}
