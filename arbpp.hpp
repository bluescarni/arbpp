#ifndef ARBPP_ARBPP_HPP
#define ARBPP_ARBPP_HPP

#include <algorithm>
#include <arb.h>
#include <arf.h>
#include <fmpr.h>
#include <iostream>
#include <limits>
#include <memory>
#include <mpfr.h>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace arbpp
{

namespace detail
{

template <typename = void>
struct base_arb
{
    // Default precision, in bits.
    static const long default_prec = 53;
};

template <typename T>
const long base_arb<T>::default_prec;

}

class arb: public detail::base_arb<>
{
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
        // Custom is_digit checker.
        static bool is_digit(char c)
        {
            const char digits[] = "0123456789";
            return std::find(digits,digits + 10,c) != (digits + 10);
        }
        // Smart pointer to handle the string output from mpfr.
        using smart_mpfr_str = std::unique_ptr<char,void (*)(char *)>;
        // Utility function to print an fmpr to stream using mpfr. Will clear f on exit.
        static void print_fmpr(std::ostream &os, ::fmpr_t f, long prec)
        {
            ::mpfr_t t;
            ::mpfr_init2(t,static_cast< ::mpfr_prec_t>(prec));
            ::fmpr_get_mpfr(t,f,MPFR_RNDN);
            // Couple of variables used below.
            const bool is_zero = (mpfr_sgn(t) == 0);
            ::mpfr_exp_t exp(0);
            char *cptr = ::mpfr_get_str(nullptr,&exp,10,0,t,MPFR_RNDN);
            // Clear everything before checking, for exception safety.
            ::mpfr_clear(t);
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
    public:
        arb() noexcept : m_prec(detail::base_arb<>::default_prec)
        {
            ::arb_init(&m_arb);
        }
        arb(const arb &other) noexcept : m_prec(other.m_prec)
        {
            ::arb_init(&m_arb);
            ::arb_set(&m_arb,other);
        }
        arb(arb &&other) noexcept : m_prec(other.m_prec)
        {
            ::arb_init(&m_arb);
            ::arb_swap(&m_arb,&other.m_arb);
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        explicit arb(T n) noexcept : arb()
        {
            ::arb_set_si(&m_arb,static_cast<long>(n));
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        explicit arb(T n) noexcept : arb()
        {
            ::arb_set_ui(&m_arb,static_cast<unsigned long>(n));
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        explicit arb(T x) noexcept : arb()
        {
            ::arf_t tmp_arf;
            ::arf_init_set_ui(tmp_arf,0u);
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_set_arf(&m_arb,tmp_arf);
            ::arf_clear(tmp_arf);
        }
        void add_error(double err) noexcept
        {
            ::arf_t tmp_arf;
            ::arf_init_set_ui(tmp_arf,0u);
            ::arf_set_d(tmp_arf,err);
            ::arb_add_error_arf(&m_arb,tmp_arf);
            ::arf_clear(tmp_arf);
        }
        ~arb()
        {
            ::arb_clear(&m_arb);
        }
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
        friend arb operator+(const arb &a, const arb &b)
        {
            const long prec = std::max<long>(a.m_prec,b.m_prec);
            arb retval;
            retval.m_prec = prec;
            ::arb_add(&retval.m_arb,a,b,prec);
            return retval;
        }
        arb cos() const
        {
            arb retval;
            ::arb_cos(&retval.m_arb,*this,m_prec);
            return retval;
        }
        operator const ::arb_struct *() const noexcept
        {
            return &m_arb;
        }
        friend std::ostream &operator<<(std::ostream &os, const arb &a)
        {
            os << '[';
            // First print the arf.
            print_arf(os,arb_midref(&a.m_arb),a.m_prec);
            os << " +/- ";
            // Now print the mag.
            print_mag(os,arb_radref(&a.m_arb),a.m_prec);
            os << ']';
            return os;
        }
        friend arb cos(const arb &a)
        {
            return a.cos();
        }
    private:
        ::arb_struct    m_arb;
        long            m_prec;
};

}

#endif
