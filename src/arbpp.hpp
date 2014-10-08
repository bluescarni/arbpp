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
 *   along with Arbpp.  If not, see <http://www.gnu.org/licenses/>.        *
 ***************************************************************************/

#ifndef ARBPP_ARBPP_HPP
#define ARBPP_ARBPP_HPP

#include <algorithm>
#include <arb.h>
#include <arf.h>
#include <cassert>
#include <cmath>
#include <fmpr.h>
#include <iostream>
#include <limits>
#include <mag.h>
#include <memory>
#include <mpfr.h>
#include <stdexcept>
#include <string>
#include <type_traits>

/// Root Arbpp namespace.
namespace arbpp
{

// Namespace for implementation details.
namespace detail
{

template <typename = void>
struct base_arb
{
    // Default precision, in bits.
    static const long default_prec = 53;
    static_assert(default_prec >= MPFR_PREC_MIN && default_prec <= MPFR_PREC_MAX,
        "Invalid default precision.");
};

template <typename T>
const long base_arb<T>::default_prec;

// Basic RAII holders.
struct arf_raii
{
    arf_raii()
    {
        // This sets to zero.
        ::arf_init(&m_arf);
    }
    arf_raii(const arf_raii &) = delete;
    arf_raii(arf_raii &&) = delete;
    arf_raii &operator=(const arf_raii &) = delete;
    arf_raii &operator=(arf_raii &&) = delete;
    operator ::arf_struct *()
    {
        return &m_arf;
    }
    operator ::arf_struct const *() const
    {
        return &m_arf;
    }
    ~arf_raii()
    {
        ::arf_clear(&m_arf);
    }
    ::arf_struct m_arf;
};

struct fmpr_raii
{
    fmpr_raii()
    {
        // Sets to zero as well.
        ::fmpr_init(&m_fmpr);
    }
    fmpr_raii(const fmpr_raii &) = delete;
    fmpr_raii(fmpr_raii &&) = delete;
    fmpr_raii &operator=(const fmpr_raii &) = delete;
    fmpr_raii &operator=(fmpr_raii &&) = delete;
    operator ::fmpr_struct *()
    {
        return &m_fmpr;
    }
    operator ::fmpr_struct const *() const
    {
        return &m_fmpr;
    }
    ~fmpr_raii()
    {
        ::fmpr_clear(&m_fmpr);
    }
    ::fmpr_struct m_fmpr;
};

struct mpfr_raii
{
    explicit mpfr_raii(::mpfr_prec_t prec)
    {
        // This will set to NaN.
        ::mpfr_init2(m_mpfr,prec);
    }
    mpfr_raii(const mpfr_raii &) = delete;
    mpfr_raii(mpfr_raii &&) = delete;
    mpfr_raii &operator=(const mpfr_raii &) = delete;
    mpfr_raii &operator=(mpfr_raii &&) = delete;
    operator std::remove_extent< ::mpfr_t>::type *()
    {
        return m_mpfr;
    }
    operator std::remove_extent< ::mpfr_t>::type const *() const
    {
        return m_mpfr;
    }
    ~mpfr_raii()
    {
        ::mpfr_clear(m_mpfr);
    }
    ::mpfr_t m_mpfr;
};

}

/// Real number represented as a floating-point ball.
/**
 * \section interop Interoperability with fundamental types
 * 
 * Interoperability with the following types is provided:
 * - \p char, <tt>signed char</tt>, \p short, \p int, \p long and unsigned counterparts,
 * - \p float and \p double.
 * 
 * \section exception_safety Exception safety guarantee
 * 
 * This class provides the strong exception safety guarantee for all operations.
 * In case of memory allocation errors by
 * lower-level libraries (e.g., GMP), the program will terminate.
 * 
 * \section move_semantics Move semantics
 *
 * Move construction and move assignment will leave the moved-from object in an
 * unspecified but valid state.
 */
class arb: private detail::base_arb<>
{
        // Import locally the RAII names.
        typedef detail::arf_raii arf_raii;
        typedef detail::fmpr_raii fmpr_raii;
        typedef detail::mpfr_raii mpfr_raii;
        // Signed integers with which arb can interoperate.
        template <typename T>
        struct is_arb_int
        {
            static const bool value = std::is_same<T,signed char>::value || std::is_same<T,short>::value ||
                std::is_same<T,int>::value || std::is_same<T,long>::value ||
                (std::is_signed<char>::value && std::is_same<T,char>::value);
        };
        // Unsigned integers with which arb can interoperate.
        template <typename T>
        struct is_arb_uint
        {
            static const bool value = std::is_same<T,unsigned char>::value || std::is_same<T,unsigned short>::value ||
                std::is_same<T,unsigned>::value || std::is_same<T,unsigned long>::value ||
                (std::is_unsigned<char>::value && std::is_same<T,char>::value);
        };
        // Floating-point types with which arb can interoperate.
        template <typename T>
        struct is_arb_float
        {
            static const bool value = std::is_same<T,float>::value || std::is_same<T,double>::value;
        };
        // Interoperable types.
        template <typename T>
        struct is_interoperable
        {
            static const bool value = is_arb_int<T>::value || is_arb_uint<T>::value || is_arb_float<T>::value;
        };
        // Custom is_digit checker.
        static bool is_digit(char c)
        {
            const char digits[] = "0123456789";
            return std::find(digits,digits + 10,c) != (digits + 10);
        }
        // Smart pointer to handle the string output from mpfr.
        typedef std::unique_ptr<char,void (*)(char *)> smart_mpfr_str;
        // Utility function to print an fmpr to stream using mpfr. Will clear f on exit.
        static void print_fmpr(std::ostream &os, ::fmpr_t f, long prec)
        {
            mpfr_raii t(static_cast< ::mpfr_prec_t>(prec));
            ::fmpr_get_mpfr(t,f,MPFR_RNDN);
            // Couple of variables used below.
            const bool is_zero = (mpfr_sgn(t.m_mpfr) == 0);
            ::mpfr_exp_t exp(0);
            char *cptr = ::mpfr_get_str(nullptr,&exp,10,0,t,MPFR_RNDN);
            // Clear before checking, for exception safety.
            ::fmpr_clear(f);
            if (!cptr) {
                throw std::invalid_argument("error while converting arb to string");
            }
            smart_mpfr_str str(cptr,::mpfr_free_str);
            // Copy into C++ string.
            std::string cpp_str(str.get());
            // Insert the radix point.
            auto it = std::find_if(cpp_str.begin(),cpp_str.end(),[](char c) {return is_digit(c);});
            if (it != cpp_str.end()) {
                ++it;
                cpp_str.insert(it,'.');
                if (exp == std::numeric_limits< ::mpfr_exp_t>::min()) {
                    throw std::invalid_argument("error while converting arb to string");
                }
                --exp;
                if (exp != ::mpfr_exp_t(0) && !is_zero) {
                    cpp_str.append(std::string("e") + std::to_string(exp));
                }
            }
            os << cpp_str;
        }
        static void print_arf(std::ostream &os, const ::arf_t a, long prec)
        {
            // Go through the conversion chain arf -> fmpr -> mpfr.
            ::fmpr_t f;
            ::fmpr_init(f);
            ::arf_get_fmpr(f,a);
            print_fmpr(os,f,prec);
        }
        static void print_mag(std::ostream &os, const ::mag_t m, long prec)
        {
            // Go through the conversion chain mag -> fmpr -> mpfr.
            ::fmpr_t f;
            ::fmpr_init(f);
            ::mag_get_fmpr(f,m);
            print_fmpr(os,f,prec);
        }
        // Generic constructor.
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        void construct(const T &n)
        {
            ::arb_set_si(&m_arb,static_cast<long>(n));
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        void construct(const T &n)
        {
            ::arb_set_ui(&m_arb,static_cast<unsigned long>(n));
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        void construct(const T &x)
        {
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_set_arf(&m_arb,tmp_arf);
        }
        // Addition.
        arb &in_place_add(const arb &other)
        {
            // Work with max precision.
            if (other.m_prec > m_prec) {
                m_prec = other.m_prec;
            }
            ::arb_add(&m_arb,&m_arb,&other.m_arb,m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        arb &in_place_add(const T &n)
        {
            ::arb_add_si(&m_arb,&m_arb,static_cast<long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        arb &in_place_add(const T &n)
        {
            ::arb_add_ui(&m_arb,&m_arb,static_cast<unsigned long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        arb &in_place_add(const T &x)
        {
            arf_raii tmp_arf;
            // Set tmp_arf *exactly* to x.
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            // Add with precision m_prec.
            ::arb_add_arf(&m_arb,&m_arb,tmp_arf,m_prec);
            return *this;
        }
        static arb binary_add(const arb &a, const arb &b)
        {
            arb retval;
            // Set max precision.
            if (a.m_prec > b.m_prec) {
                retval.m_prec = a.m_prec;
            } else {
                retval.m_prec = b.m_prec;
            }
            ::arb_add(&retval.m_arb,&a.m_arb,&b.m_arb,retval.m_prec);
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_add(const arb &a, const T &n)
        {
            arb retval;
            ::arb_add_si(&retval.m_arb,&a.m_arb,static_cast<long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_add(const T &n, const arb &a)
        {
            return binary_add(a,n);
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_add(const arb &a, const T &n)
        {
            arb retval;
            ::arb_add_ui(&retval.m_arb,&a.m_arb,static_cast<unsigned long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_add(const T &n, const arb &a)
        {
            return binary_add(a,n);
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_add(const arb &a, const T &x)
        {
            arb retval;
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_add_arf(&retval.m_arb,&a.m_arb,tmp_arf,a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_add(const T &x, const arb &a)
        {
            return binary_add(a,x);
        }
        // Subtraction.
        arb &in_place_sub(const arb &other)
        {
            // Work with max precision.
            if (other.m_prec > m_prec) {
                m_prec = other.m_prec;
            }
            ::arb_sub(&m_arb,&m_arb,&other.m_arb,m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        arb &in_place_sub(const T &n)
        {
            ::arb_sub_si(&m_arb,&m_arb,static_cast<long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        arb &in_place_sub(const T &n)
        {
            ::arb_sub_ui(&m_arb,&m_arb,static_cast<unsigned long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        arb &in_place_sub(const T &x)
        {
            arf_raii tmp_arf;
            // Set tmp_arf *exactly* to x.
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            // Sub with precision m_prec.
            ::arb_sub_arf(&m_arb,&m_arb,tmp_arf,m_prec);
            return *this;
        }
        static arb binary_sub(const arb &a, const arb &b)
        {
            arb retval;
            // Set max precision.
            if (a.m_prec > b.m_prec) {
                retval.m_prec = a.m_prec;
            } else {
                retval.m_prec = b.m_prec;
            }
            ::arb_sub(&retval.m_arb,&a.m_arb,&b.m_arb,retval.m_prec);
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_sub(const arb &a, const T &n)
        {
            arb retval;
            ::arb_sub_si(&retval.m_arb,&a.m_arb,static_cast<long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_sub(const T &n, const arb &a)
        {
            auto retval = binary_sub(a,n);
            retval.negate();
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_sub(const arb &a, const T &n)
        {
            arb retval;
            ::arb_sub_ui(&retval.m_arb,&a.m_arb,static_cast<unsigned long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_sub(const T &n, const arb &a)
        {
            auto retval = binary_sub(a,n);
            retval.negate();
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_sub(const arb &a, const T &x)
        {
            arb retval;
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_sub_arf(&retval.m_arb,&a.m_arb,tmp_arf,a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_sub(const T &x, const arb &a)
        {
            auto retval = binary_sub(a,x);
            retval.negate();
            return retval;
        }
        // Multiplication.
        arb &in_place_mul(const arb &other)
        {
            if (other.m_prec > m_prec) {
                m_prec = other.m_prec;
            }
            ::arb_mul(&m_arb,&m_arb,&other.m_arb,m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        arb &in_place_mul(const T &n)
        {
            ::arb_mul_si(&m_arb,&m_arb,static_cast<long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        arb &in_place_mul(const T &n)
        {
            ::arb_mul_ui(&m_arb,&m_arb,static_cast<unsigned long>(n),m_prec);
            return *this;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        arb &in_place_mul(const T &x)
        {
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_mul_arf(&m_arb,&m_arb,tmp_arf,m_prec);
            return *this;
        }
        static arb binary_mul(const arb &a, const arb &b)
        {
            arb retval;
            // Set max precision.
            if (a.m_prec > b.m_prec) {
                retval.m_prec = a.m_prec;
            } else {
                retval.m_prec = b.m_prec;
            }
            ::arb_mul(&retval.m_arb,&a.m_arb,&b.m_arb,retval.m_prec);
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_mul(const arb &a, const T &n)
        {
            arb retval;
            ::arb_mul_si(&retval.m_arb,&a.m_arb,static_cast<long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        static arb binary_mul(const T &n, const arb &a)
        {
            return binary_mul(a,n);
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_mul(const arb &a, const T &n)
        {
            arb retval;
            ::arb_mul_ui(&retval.m_arb,&a.m_arb,static_cast<unsigned long>(n),a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        static arb binary_mul(const T &n, const arb &a)
        {
            return binary_mul(a,n);
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_mul(const arb &a, const T &x)
        {
            arb retval;
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_mul_arf(&retval.m_arb,&a.m_arb,tmp_arf,a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_mul(const T &x, const arb &a)
        {
            return binary_mul(a,x);
        }
        // Enabler for the generic ctor.
        template <typename T>
        using generic_enabler = typename std::enable_if<is_interoperable<T>::value,int>::type;
    public:
        /// Default precision.
        /**
         * @return default precision for arbpp::arb objects.
         */
        static long get_default_precision()
        {
            return detail::base_arb<>::default_prec;
        }
        /// Default constructor.
        /**
         * The value is initialised with both midpoint and radius equal to zero.
         * The precision is set to the default value.
         */
        arb() : m_prec(get_default_precision())
        {
            ::arb_init(&m_arb);
        }
        /// Copy constructor.
        /**
         * @param[in] other construction argument.
         */
        arb(const arb &other) : m_prec(other.m_prec)
        {
            ::arb_init(&m_arb);
            ::arb_set(&m_arb,&other.m_arb);
        }
        /// Move constructor.
        /**
         * @param[in] other construction argument.
         */
        arb(arb &&other) noexcept : m_prec(get_default_precision())
        {
            // Init a default arb, swap it out with other.
            ::arb_init(&m_arb);
            swap(other);
        }
        /// Generic constructor.
        /**
         * \note
         * This constructor is enabled only if \p T is an \ref interop "interoperable type".
         * 
         * Construction from an interoperable type is always exact. The precision
         * will be set to the default value.
         * 
         * @param[in] x construction argument.
         */
        template <typename T, generic_enabler<T> = 0>
        explicit arb(const T &x) : m_prec(get_default_precision())
        {
            ::arb_init(&m_arb);
            construct(x);
        }
        /// Generic constructor with precision.
        /**
         * \note
         * This constructor is enabled only if \p T is an \ref interop "interoperable type".
         * 
         * This constructor will first initialise \p this to \p x, and then it will round \p this
         * to the precision \p prec. The precision of \p this will also be set to \p prec.
         * 
         * @param[in] x construction argument.
         * @param[in] prec desired precision.
         * 
         * @throws unspecified any exception thrown by arb::set_precision().
         */
        template <typename T, generic_enabler<T> = 0>
        explicit arb(const T &x, long prec)
        {
            // Set the precision value, with error checking.
            set_precision(prec);
            // Construct as usual.
            ::arb_init(&m_arb);
            construct(x);
            // Round-set self.
            ::arb_set_round(&m_arb,&m_arb,m_prec);
        }
        /// Constructor from string.
        /**
         * This constructor will set the midpoint to a value represented
         * by the string \p str and rounded to a precision of \p prec bits.
         * The expected string format is the same described in the MPFR documentation,
         * and base 10 representation is assumed.
         * 
         * If the value of the input string can be represented exactly with \p prec
         * bits of precision, then the midpoint will be set to that value
         * and the radius will be set to zero. Otherwise, the midpoint will be set
         * to the nearest representable value and the radius is guaranteed to include
         * the original value.
         * 
         * Please note that the entire string must represent a valid floating-point
         * value.
         * 
         * @param[in] str string used for construction.
         * @param[in] prec desired precision.
         * 
         * @throws std::invalid_argument in case of an invalid input string.
         * @throws std::underflow_error in case of underflow.
         * @throws unspecified any exception thrown by arb::set_precision().
         */
        arb(const std::string &str, long prec = arb::get_default_precision())
        {
            // Clear the MPFR underflow flag.
            ::mpfr_clear_underflow();
            // Try to parse an mpfr from the input string.
            mpfr_raii m(prec);
            char *endptr;
            const int retval = ::mpfr_strtofr(m,str.c_str(),&endptr,10,MPFR_RNDN);
            // NOTE: the first case corresponds to invalid string input, the second
            // to valid input only for a subset of the string.
            if (endptr == str.c_str() || endptr != str.c_str() + str.size()) {
                throw std::invalid_argument("invalid string input");
            }
            if (::mpfr_underflow_p()) {
                ::mpfr_clear_underflow();
                throw std::underflow_error("underflow in conversion from string");
            }
            // Set precision.
            set_precision(prec);
            // NOTE: here we could have non-finite number in the following cases:
            // - the passed string was already representing a non-finite ("inf", "nan", etc.),
            // - the passed string represents a number which gets rounded to +-inf.
            // In both cases it does not make sense to do the radius calculation, just leave
            // the non-finite midpoint and a zero radius.
            if (retval && ::mpfr_number_p(m)) {
                mpfr_raii tmp0(prec), tmp1(prec);
                ::mpfr_set(tmp0.m_mpfr,m.m_mpfr,MPFR_RNDN);
                ::mpfr_set(tmp1.m_mpfr,m.m_mpfr,MPFR_RNDN);
                if (retval < 0) {
                    // Converted value is lower than the exact value.
                    ::mpfr_nextabove(tmp0);
                    // NOTE: it seems like the behaviour of mpfr when
                    // the exact value is above the max or below the min is always to round to +-inf,
                    // so we should never enter here in a situation in which going above/below would lead
                    // to a nonfinite value.
                    assert(::mpfr_number_p(tmp0));
                    ::mpfr_sub(tmp1,tmp0,m,MPFR_RNDN);
                } else {
                    // Converted value is higher than the exact value.
                    ::mpfr_nextbelow(tmp0);
                    assert(::mpfr_number_p(tmp0));
                    ::mpfr_sub(tmp1,m,tmp0,MPFR_RNDN);
                }
                // Check if any underflow happened.
                if (::mpfr_underflow_p()) {
                    ::mpfr_clear_underflow();
                    throw std::underflow_error("underflow in the calculation of the radius");
                }
                // Construct an empty instance and set the midpoint.
                ::arb_init(&m_arb);
                ::arf_set_mpfr(arb_midref(&m_arb),m);
                // Set the radius via a temporary fmpr.
                fmpr_raii f;
                ::fmpr_set_mpfr(f,tmp1);
                ::mag_set_fmpr(arb_radref(&m_arb),f);
            } else {
                // Just construct empty instance and set the midpoint.
                // NOTE: we need to repeat these two lines here for exception safety (instead
                // of sharing them outside the if block).
                ::arb_init(&m_arb);
                ::arf_set_mpfr(arb_midref(&m_arb),m);
            }
        }
        /// Destructor.
        ~arb()
        {
            ::arb_clear(&m_arb);
        }
        /// Copy assignment.
        /**
         * @param[in] other assignment argument.
         * 
         * @return reference to \p this.
         */
        arb &operator=(const arb &other)
        {
            if (this == &other) {
                return *this;
            }
            ::arb_set(&m_arb,&other.m_arb);
            m_prec = other.m_prec;
            return *this;
        }
        /// Move assignment.
        /**
         * @param[in] other assignment argument.
         * 
         * @return reference to \p this.
         */
        arb &operator=(arb &&other) noexcept
        {
            // this == &other check already in swap().
            swap(other);
            return *this;
        }
        /// Generic assignment.
        /**
         * \note
         * This assignment operator is enabled only if \p T is an \ref interop "interoperable type".
         * 
         * The operation is equivalent to an assignment from an arbpp::arb object constructed from \p x.
         * 
         * @param[in] x assignment argument.
         * 
         * @return reference to \p this.
         */
        template <typename T, generic_enabler<T> = 0>
        arb &operator=(const T &x)
        {
            m_prec = get_default_precision();
            construct(x);
            return *this;
        }
        /// Add error.
        /**
         * Add \p err to the radius.
         * 
         * @param[in] err error value.
         * 
         * @throws std::invalid_argument if \p err is negative or NaN.
         */
        void add_error(double err)
        {
            if (std::isnan(err) || err <= 0.) {
                throw std::invalid_argument("an error value must be positive and not NaN");
            }
            arf_raii tmp_arf;
            ::arf_set_d(tmp_arf,err);
            ::arb_add_error_arf(&m_arb,tmp_arf);
        }
        /// Precision setter.
        /**
         * Set the precision of \p this to \p prec bits.
         * 
         * @param[in] prec desired value for the precision.
         * 
         * @throws std::invalid_argument if \p prec is not positive or not within
         * an implementation-defined range.
         */
        void set_precision(long prec)
        {
            // NOTE: this precision is to be used within mpfr routines as well, so check that
            // the value we are setting is not outside the mpfr bounds.
            // The mpfr headers say that mpfr_prec_t is always a signed int,
            // so there is no danger here in comparing to long.
            if (prec < 1 || prec < MPFR_PREC_MIN || prec > MPFR_PREC_MAX) {
                throw std::invalid_argument("invalid precision value");
            }
            m_prec = prec;
        }
        /// Precision getter.
        /**
         * @return precision associated to \p this.
         */
        long get_precision() const
        {
            return m_prec;
        }
        /// Swap method.
        /**
         * Swap efficiently \p this with \p other.
         * 
         * @param[in] other argument for swap.
         */
        void swap(arb &other) noexcept
        {
            if (this == &other) {
                return;
            }
            ::arb_swap(&m_arb,&other.m_arb);
            std::swap(m_prec,other.m_prec);
        }
        /// Get a const pointer to the internal \p arb_struct.
        /**
         * @return const pointer to the internal \p arb_struct.
         */
        const ::arb_struct *get_arb_t() const
        {
            return &m_arb;
        }
        /// Get a mutable pointer to the internal \p arb_struct.
        /**
         * @return pointer to the internal \p arb_struct.
         */
        ::arb_struct *get_arb_t()
        {
            return &m_arb;
        }
        /// Stream operator.
        /**
         * This function will print to stream a human-readable representation
         * of \p a.
         * 
         * @param[in,out] os target stream.
         * @param[in] a arbpp::arb to be streamed.
         * 
         * @return reference to \p os.
         * 
         * @throws std::invalid_argument in case of any error in the conversion
         * of \p a to string.
         */
        friend std::ostream &operator<<(std::ostream &os, const arb &a)
        {
            os << '(';
            // First print the arf.
            print_arf(os,arb_midref(&a.m_arb),a.m_prec);
            os << " +/- ";
            // Now print the mag.
            print_mag(os,arb_radref(&a.m_arb),a.m_prec);
            os << ')';
            return os;
        }
        /// Midpoint getter.
        /**
         * @return the midpoint of \p this, rounded to \p double in an unspecified
         * direction.
         */
        double get_midpoint() const
        {
            return ::arf_get_d(arb_midref(&m_arb),ARF_RND_DOWN);
        }
        /// Radius getter.
        /**
         * @return the radius of \p this, rounded to \p double in an unspecified
         * direction.
         */
        double get_radius() const
        {
            fmpr_raii tmp;
            ::mag_get_fmpr(tmp,arb_radref(&m_arb));
            return ::fmpr_get_d(tmp,FMPR_RND_UP);
        }
        /// Identity operator.
        /**
         * @return a copy of \p this.
         */
        arb operator+() const
        {
            return *this;
        }
        /// In-place addition.
        /**
         * \note
         * This operator is enabled only if \p T is an \ref interop "interoperable type"
         * or arbpp::arb.
         * 
         * This method will set \p this to <tt>this + x</tt>. In case \p T is arbpp::arb, then
         * the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p this and \p x. Otherwise, the operation is carried out with the
         * precision of \p this.
         * 
         * @param[in] x addition argument.
         * 
         * @return reference to \p this.
         */
        template <typename T>
        auto operator+=(const T &x) -> decltype(this->in_place_add(x))
        {
            return in_place_add(x);
        }
        /// Generic binary addition involving arbpp::arb.
        /**
         * \note
         * This template operator is enabled only if either:
         * - \p T is arbpp::arb and \p U is an \ref interop "interoperable type",
         * - \p U is arbpp::arb and \p T is an \ref interop "interoperable type",
         * - both \p T and \p U are arbpp::arb.
         * 
         * This method will compute <tt>a + b</tt> and return it as an arbpp::arb. In case \p T and \p U are both
         * arbpp::arb, then the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p a and \p b. Otherwise, the result will have the precision of the
         * arbpp::arb argument.
         * 
         * @param[in] a first operand.
         * @param[in] b second operand.
         * 
         * @return <tt>a + b</tt>.
         */
        template <typename T, typename U>
        friend auto operator+(const T &a, const U &b) -> decltype(arb::binary_add(a,b))
        {
            return binary_add(a,b);
        }
        /// Negation.
        /**
         * Negate in-place.
         */
        void negate()
        {
            ::arb_neg(&m_arb,&m_arb);
        }
        /// Negated copy.
        /**
         * @return a negated copy of \p this.
         */
        arb operator-() const
        {
            arb retval{*this};
            retval.negate();
            return retval;
        }
        /// In-place subtraction.
        /**
         * \note
         * This operator is enabled only if \p T is an \ref interop "interoperable type"
         * or arbpp::arb.
         * 
         * This method will set \p this to <tt>this - x</tt>. In case \p T is arbpp::arb, then
         * the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p this and \p x. Otherwise, the operation is carried out with the
         * precision of \p this.
         * 
         * @param[in] x subtraction argument.
         * 
         * @return reference to \p this.
         */
        template <typename T>
        auto operator-=(const T &x) -> decltype(this->in_place_sub(x))
        {
            return in_place_sub(x);
        }
        /// Generic binary subtraction involving arbpp::arb.
        /**
         * \note
         * This template operator is enabled only if either:
         * - \p T is arbpp::arb and \p U is an \ref interop "interoperable type",
         * - \p U is arbpp::arb and \p T is an \ref interop "interoperable type",
         * - both \p T and \p U are arbpp::arb.
         * 
         * This method will compute <tt>a - b</tt> and return it as an arbpp::arb. In case \p T and \p U are both
         * arbpp::arb, then the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p a and \p b. Otherwise, the result will have the precision of the
         * arbpp::arb argument.
         * 
         * @param[in] a first operand.
         * @param[in] b second operand.
         * 
         * @return <tt>a - b</tt>.
         */
        template <typename T, typename U>
        friend auto operator-(const T &a, const U &b) -> decltype(arb::binary_sub(a,b))
        {
            return binary_sub(a,b);
        }
        /// In-place multiplication.
        /**
         * \note
         * This operator is enabled only if \p T is an \ref interop "interoperable type"
         * or arbpp::arb.
         * 
         * This method will set \p this to <tt>this * x</tt>. In case \p T is arbpp::arb, then
         * the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p this and \p x. Otherwise, the operation is carried out with the
         * precision of \p this.
         * 
         * @param[in] x multiplication argument.
         * 
         * @return reference to \p this.
         */
        template <typename T>
        auto operator*=(const T &x) -> decltype(this->in_place_mul(x))
        {
            return in_place_mul(x);
        }
        /// Generic binary multiplication involving arbpp::arb.
        /**
         * \note
         * This template operator is enabled only if either:
         * - \p T is arbpp::arb and \p U is an \ref interop "interoperable type",
         * - \p U is arbpp::arb and \p T is an \ref interop "interoperable type",
         * - both \p T and \p U are arbpp::arb.
         * 
         * This method will compute <tt>a * b</tt> and return it as an arbpp::arb. In case \p T and \p U are both
         * arbpp::arb, then the operation is carried out with a precision corresponding to the maximum between
         * the precisions of \p a and \p b. Otherwise, the result will have the precision of the
         * arbpp::arb argument.
         * 
         * @param[in] a first operand.
         * @param[in] b second operand.
         * 
         * @return <tt>a * b</tt>.
         */
        template <typename T, typename U>
        friend auto operator*(const T &a, const U &b) -> decltype(arb::binary_mul(a,b))
        {
            return binary_mul(a,b);
        }
        /// Cosine.
        /**
         * @return the cosine of \p this.
         */
        arb cos() const
        {
            arb retval;
            ::arb_cos(&retval.m_arb,&m_arb,m_prec);
            retval.m_prec = m_prec;
            return retval;
        }
    private:
        ::arb_struct    m_arb;
        long            m_prec;
};

/// Cosine.
/**
 * @param[in] a cosine argument.
 * 
 * @return <tt>a.cos()</tt>.
 */
inline arb cos(const arb &a)
{
    return a.cos();
}

/// Swap.
/**
 * Equivalent to <tt>a0.swap(a1)</tt>.
 * 
 * @param[in] a0 first argument.
 * @param[in] a1 second argument.
 */
inline void swap(arb &a0, arb &a1) noexcept
{
    a0.swap(a1);
}

/// Literal namespace.
inline namespace literals
{

/// Literal for arbpp::arb.
/**
 * The return value will be constructed with the default precision.
 *
 * @param[in] s literal string.
 *
 * @return an arbpp::arb constructed from \p s.
 *
 * @throws unspecified any exception thrown by the constructor of
 * arbpp::arb from string.
 */
inline arb operator "" _arb(const char *s)
{
    return arb(s);
}

}

}

#endif
