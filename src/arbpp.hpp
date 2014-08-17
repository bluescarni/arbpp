#ifndef ARBPP_ARBPP_HPP
#define ARBPP_ARBPP_HPP

#include <algorithm>
#include <arb.h>
#include <arf.h>
#include <cmath>
#include <fmpr.h>
#include <iostream>
#include <limits>
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
};

template <typename T>
const long base_arb<T>::default_prec;

}

// A few forward declarations.
class arb;
arb operator+(const arb &, const arb &);
std::ostream &operator<<(std::ostream &, const arb &);

/// Real number represented as a floating-point ball.
/**
 * \section interop Interoperability with fundamental types
 * 
 * Interoperability with the following types is provided:
 * - <tt>signed char</tt>, \p short, \p int, \p long and unsigned counterparts,
 * - \p char,
 * - \p float and \p double.
 * 
 * \section move_semantics Move semantics
 *
 * Move construction and move assignment will leave the moved-from object in an
 * unspecified but valid state.
 */
class arb: public detail::base_arb<>
{
        // Friends.
        friend arb operator+(const arb &, const arb &);
        friend std::ostream &operator<<(std::ostream &, const arb &);
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
            ::arf_t tmp_arf;
            ::arf_init_set_ui(tmp_arf,0u);
            ::arf_set_d(tmp_arf,static_cast<double>(x));
            ::arb_set_arf(&m_arb,tmp_arf);
            ::arf_clear(tmp_arf);
        }
    public:
        /// Default constructor.
        /**
         * The value is initialised with both midpoint and radius equal to zero.
         */
        arb() noexcept : m_prec(detail::base_arb<>::default_prec)
        {
            ::arb_init(&m_arb);
        }
        /// Copy constructor.
        /**
         * @param[in] other construction argument.
         */
        arb(const arb &other) noexcept : m_prec(other.m_prec)
        {
            ::arb_init(&m_arb);
            ::arb_set(&m_arb,other);
        }
        /// Move constructor.
        /**
         * @param[in] other construction argument.
         */
        arb(arb &&other) noexcept : m_prec(other.m_prec)
        {
            ::arb_init(&m_arb);
            ::arb_swap(&m_arb,&other.m_arb);
        }
        /// Generic constructor.
        /**
         * \note
         * This constructor is enabled only if \p T is an \ref interop "interoperable type".
         * 
         * Construction from an interoperable type is always exact. The internal precision
         * of the object will be set to the default value.
         * 
         * @param[in] x construction argument.
         */
        template <typename T, typename std::enable_if<is_interoperable<T>::value,int>::type = 0>
        explicit arb(T x) noexcept : m_prec(detail::base_arb<>::default_prec)
        {
            ::arb_init(&m_arb);
            construct(x);
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
        arb &operator=(const arb &other) noexcept
        {
            if (this == &other) {
                return *this;
            }
            ::arb_set(&m_arb,other);
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
            swap(other);
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
            if (err < 0. || std::isnan(err)) {
                throw std::invalid_argument("an error value must be non-negative and not NaN");
            }
            ::arf_t tmp_arf;
            ::arf_init_set_ui(tmp_arf,0u);
            ::arf_set_d(tmp_arf,err);
            ::arb_add_error_arf(&m_arb,tmp_arf);
            ::arf_clear(tmp_arf);
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
        long get_precision() const noexcept
        {
            return m_prec;
        }
        /// Cosine.
        /**
         * @return the cosine of \p this.
         */
        arb cos() const noexcept
        {
            arb retval;
            retval.m_prec = m_prec;
            ::arb_cos(&retval.m_arb,*this,m_prec);
            return retval;
        }
        /// Swap method.
        void swap(arb &other) noexcept
        {
            if (this == &other) {
                return;
            }
            ::arb_swap(&m_arb,&other.m_arb);
            std::swap(m_prec,other.m_prec);
        }
        /// Implicit conversion operator to <tt>const arb_struct *</tt>.
        /**
         * This operator allows to pass an arbpp::arb object
         * as a <tt>const arb_t</tt> argument to Arb's C API.
         */
        operator const ::arb_struct *() const noexcept
        {
            return &m_arb;
        }
    private:
        ::arb_struct    m_arb;
        long            m_prec;
};

/// Binary addition.
/**
 * @param[in] a first operand.
 * @param[in] b second operand.
 * 
 * @return <tt>a+b</tt>.
 */
inline arb operator+(const arb &a, const arb &b)
{
    const long prec = std::max<long>(a.m_prec,b.m_prec);
    arb retval;
    retval.m_prec = prec;
    ::arb_add(&retval.m_arb,a,b,prec);
    return retval;
}

/// Stream operator.
inline std::ostream &operator<<(std::ostream &os, const arb &a)
{
    os << '[';
    // First print the arf.
    arb::print_arf(os,arb_midref(&a.m_arb),a.m_prec);
    os << " +/- ";
    // Now print the mag.
    arb::print_mag(os,arb_radref(&a.m_arb),a.m_prec);
    os << ']';
    return os;
}

/// Cosine.
inline arb cos(const arb &a)
{
    return a.cos();
}

/// Swap.
inline void swap(arb &a0, arb &a1)
{
    a0.swap(a1);
}

}

#endif
