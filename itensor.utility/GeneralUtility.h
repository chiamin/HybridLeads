#ifndef __GENERALUTILITY_H_CMC__
#define __GENERALUTILITY_H_CMC__

// Use iprint
#ifndef iprint
#define iprint(name) myprinter(#name, (name))
#endif

using namespace std;

template <typename TType>
inline void myprinter (string name, const TType& value)
{
    cout << name << endl;
    cout << value << endl;
}


// Variadic function Template that takes
// variable number of arguments and prints
// all of them.
void variadic_print () {} // Terminate function
template <typename T, typename... Types>
void variadic_print(T var1, Types... var2)
{
    cout << var1 << " ";
    variadic_print (var2...);
}

#define mycheck(condition, message) mycheck_impl(condition, __func__, message)
template <typename TypeName>
inline void mycheck_impl (const TypeName& condition, const string& func_name, string message="")
{
    if (!bool(condition))
    {
        cout << func_name << ": " << message << endl;
        throw;
    }
}
/*
template <typename ConditionT, typename... MessageT>
inline void mycheck (const ConditionT& condition, const MessageT& ...message,
                     const std::source_location& location = std::source_location::current())
{
    if (!bool(condition))
    {
        cout << location.function_name() << ": ";
        variadic_print (message...);
        cout << endl;
        throw;
    }
}*/

namespace iutility {

inline complex<double> conjT (complex<double> a)
{
    return std::conj(a);
}

inline double conjT (double a)
{
    return a;
}

}

namespace iut {

inline complex<double> conj (complex<double> a)
{
    return std::conj(a);
}

inline double conj (double a)
{
    return a;
}

}
#endif
