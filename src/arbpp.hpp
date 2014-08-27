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

// Basic RAII holder for arf objects.
struct arf_raii
{
    arf_raii()
    {
        ::arf_init_set_ui(&m_arf,0u);
    }
    arf_raii(const arf_raii &) = delete;
    arf_raii(arf_raii &&) = delete;
    arf_raii &operator=(const arf_raii &) = delete;
    arf_raii &operator=(arf_raii &&) = delete;
    ~arf_raii()
    {
        ::arf_clear(&m_arf);
    }
    ::arf_struct m_arf;
};

}

/// Real number represented as a floating-point ball.
/**
 * \section interop Interoperability with fundamental types
 * 
 * Interoperability with the following types is provided:
 * - <tt>signed char</tt>, \p short, \p int, \p long and unsigned counterparts,
 * - \p char, \p float and \p double.
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
class arb: public detail::base_arb<>
{
        typedef detail::arf_raii arf_raii;
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
        // Type valid for unary arithmetic operations.
        template <typename T>
        struct is_unary_arithmetic_type
        {
            static const bool value = is_interoperable<T>::value || std::is_same<T,arb>::value;
        };
        // Type valid for binary arithmetic operations.
        template <typename T, typename U>
        struct are_binary_arithmetic_types
        {
            static const bool value = (std::is_same<T,arb>::value && is_interoperable<U>::value) ||
                (std::is_same<U,arb>::value && is_interoperable<T>::value) ||
                (std::is_same<T,arb>::value && std::is_same<U,arb>::value);
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
            arf_raii tmp_arf;
            ::arf_set_d(&tmp_arf.m_arf,static_cast<double>(x));
            ::arb_set_arf(&m_arb,&tmp_arf.m_arf);
        }
        // Addition.
        void in_place_add(const arb &other)
        {
            // Work with max precision.
            if (other.m_prec > m_prec) {
                m_prec = other.m_prec;
            }
            ::arb_add(&m_arb,&m_arb,&other.m_arb,m_prec);
        }
        template <typename T, typename std::enable_if<is_arb_int<T>::value,int>::type = 0>
        void in_place_add(const T &n)
        {
            ::arb_add_si(&m_arb,&m_arb,static_cast<long>(n),m_prec);
        }
        template <typename T, typename std::enable_if<is_arb_uint<T>::value,int>::type = 0>
        void in_place_add(const T &n)
        {
            ::arb_add_ui(&m_arb,&m_arb,static_cast<unsigned long>(n),m_prec);
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        void in_place_add(const T &x)
        {
            arf_raii tmp_arf;
            ::arf_set_d(&tmp_arf.m_arf,static_cast<double>(x));
            ::arb_add_arf(&m_arb,&m_arb,&tmp_arf.m_arf,m_prec);
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
            ::arf_set_d(&tmp_arf.m_arf,static_cast<double>(x));
            ::arb_add_arf(&retval.m_arb,&a.m_arb,&tmp_arf.m_arf,a.m_prec);
            retval.m_prec = a.m_prec;
            return retval;
        }
        template <typename T, typename std::enable_if<is_arb_float<T>::value,int>::type = 0>
        static arb binary_add(const T &x, const arb &a)
        {
            return binary_add(a,x);
        }
    public:
        /// Default constructor.
        /**
         * The value is initialised with both midpoint and radius equal to zero.
         * The precision is set to the default value.
         */
        arb() : m_prec(detail::base_arb<>::default_prec)
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
        arb(arb &&other) noexcept : m_prec(detail::base_arb<>::default_prec)
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
        template <typename T, typename std::enable_if<is_interoperable<T>::value,int>::type = 0>
        explicit arb(T x) : m_prec(detail::base_arb<>::default_prec)
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
        arb &operator=(arb &&other)
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
            if (std::isnan(err) || err <= 0.) {
                throw std::invalid_argument("an error value must be positive and not NaN");
            }
            arf_raii tmp_arf;
            ::arf_set_d(&tmp_arf.m_arf,err);
            ::arb_add_error_arf(&m_arb,&tmp_arf.m_arf);
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
         * The number of decimal digits that will be displayed depends on
         * the precision value associated to \p a.
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
         * the precisions of \p this and \p x.
         * 
         * @param[in] x addition argument.
         * 
         * @return reference to \p this.
         */
        template <typename T, typename std::enable_if<is_unary_arithmetic_type<T>::value,int>::type = 0>
        arb &operator+=(const T &x)
        {
            in_place_add(x);
            return *this;
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
         * @return <tt>a+b</tt>.
         */
        template <typename T, typename U, typename std::enable_if<are_binary_arithmetic_types<T,U>::value,int>::type = 0>
        friend arb operator+(const T &a, const U &b)
        {
            return binary_add(a,b);
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

}

#endif
