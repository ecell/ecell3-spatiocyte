#ifndef PACKER_HPP
#define PACKER_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_integral.hpp>

template<typename T>
struct storage_for
{
    struct type
    {
        type(T const& data)
        {
            *this = *reinterpret_cast<T const*>(&data);
        }

        unsigned char _[sizeof(T)];
    };
};

struct be_packer
{
    template<typename T>
    std::size_t operator()(unsigned char *buf, T const& data, typename boost::enable_if<boost::is_integral<T> >::type* = 0)
    {
        pack_integral<T, 0>(buf, data);
        return sizeof(T);
    }

    template<typename T>
    std::size_t operator()(unsigned char *buf, T const& data, typename boost::disable_if<boost::is_integral<T> >::type* = 0)
    {
        pack_float<T, 0>(buf, data);
        return sizeof(T);
    }

private:
    template<typename T, unsigned P>
    static typename boost::enable_if<boost::mpl::bool_<P == sizeof(T) * 8> >::type* pack_integral(unsigned char*, T const&)
    {
        return 0;
    }

    template<typename T, unsigned P>
    static typename boost::disable_if<boost::mpl::bool_<P == sizeof(T) * 8> >::type* pack_integral(unsigned char* p, T const& data)
    {
        *p = ((data >> ((sizeof(T) - 1) * 8 - P)) & 0xff);
        pack_integral<T, P + 8>(p + 1, data);
        return 0;
    }

    template<typename T, unsigned P>
    static void pack_float(unsigned char* p, T const& data)
    {
        typedef typename storage_for<T>::type storage;
        static const int a = 0x01020304;
        if (reinterpret_cast<unsigned char const*>(&a)[0] == 0x01)
        {
            *reinterpret_cast<storage*>(p) = *reinterpret_cast<storage const*>(&data);
        }
        else
        {
            for (std::size_t i = 0; i < sizeof(T); ++i)
            {
                p[i] = reinterpret_cast<unsigned char*>(p)[(sizeof(T) - 1) - i];
            }
        }
    }
};

struct le_packer
{
    template<typename T>
    std::size_t operator()(unsigned char *buf, T const& data, typename boost::enable_if<boost::is_integral<T> >::type* = 0)
    {
        pack_integral<T, 0>(buf, data);
        return sizeof(T);
    }

    template<typename T>
    std::size_t operator()(unsigned char *buf, T const& data, typename boost::disable_if<boost::is_integral<T> >::type* = 0)
    {
        pack_float<T, 0>(buf, data);
        return sizeof(T);
    }

private:
    template<typename T, unsigned P>
    static typename boost::enable_if<boost::mpl::bool_<P == sizeof(T) * 8> >::type* pack_integral(unsigned char*, T const&)
    {
        return 0;
    }

    template<typename T, unsigned P>
    static typename boost::disable_if<boost::mpl::bool_<P == sizeof(T) * 8> >::type* pack_integral(unsigned char* p, T const& data)
    {
        *p = ((data >> P) & 0xff);
        pack_integral<T, P + 8>(p + 1, data);
        return 0;
    }

    template<typename T, unsigned P>
    static void pack_float(unsigned char* p, T const& data)
    {
        typedef typename storage_for<T>::type storage;
        static const int a = 0x01020304;
        if (reinterpret_cast<unsigned char const*>(&a)[0] == 0x04)
        {
            *reinterpret_cast<storage*>(p) = *reinterpret_cast<storage const*>(&data);
        }
        else
        {
            for (std::size_t i = 0; i < sizeof(T); ++i)
            {
                p[i] = reinterpret_cast<unsigned char*>(p)[(sizeof(T) - 1) - i];
            }
        }
    }
};

#endif /* PACKER_HPP */
