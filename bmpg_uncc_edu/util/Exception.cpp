////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_exception_cpp
#define bmpg_uncc_edu_util_exception_cpp

#include <sstream>
#include <iostream>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace util {

Exception::Exception() : message_(""), file_(NULL), line_(-1)
{
}

Exception::Exception(const char * msg, const char * file, int line) : 
	message_(msg), file_(file), line_(line)
{
}

Exception::Exception(const std::string & msg, const char * file, int line) :
	message_(msg), file_(file), line_(line)
{
}

Exception::~Exception() throw ()
{
}
		
inline const char * Exception::what() const throw()
{
	const char * msg = "Exception occured - default message";
	try {
		std::stringstream sstr;
		sstr << "Exception";// type = " << typeid((*this)).name();
		if (file_ != NULL) {
			sstr << " : In file \"" << file_ << "\"";
			if (line_ >= 0)
				sstr << " at line number " << line_;
		}
		if (!message_.empty())
			sstr << " : message = \"" << message_.c_str() << "\"";	

		msg = sstr.str().c_str();
	}
	catch (exception & e) {
		// ignore any exception that might be thrown in the process of 
		// formatting the message.  In this case, the default message will
		// be retunred.
	}
	return msg;
}

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif

