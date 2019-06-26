////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_exception_hpp
#define bmpg_uncc_edu_util_exception_hpp

#include <exception>
#include <string>

namespace bmpg_uncc_edu {
namespace util{

class Exception : public std::exception
{
	public:
		Exception();
		Exception(const char * message, const char * file = NULL, int line = -1);
		Exception(const std::string &, const char * file = NULL, int line = -1);
		~Exception() throw ();

		const char * what() const throw();
	
	private:
		std::string message_;
		const char * file_;
		int line_;
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif

