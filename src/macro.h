#ifndef CWDT_MACRO_H
#define CWDT_MACRO_H

#ifdef CWDT_VERBOSE
#	define VERBOSE(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << x << std::endl // [yellow] + white cout
#	define WARNING(x) std::cerr << "\033[33m" << "[" << __FILE__ << " " << __LINE__ << "]" << "\033[0m" << " " << "\033[31m" << x << "\033[0m" << std::endl // [yellow] + red cout
#else
#	define VERBOSE(x)
#	define WARNING(x)
#endif

#ifdef CWDT_DEBUG
#endif

#endif
