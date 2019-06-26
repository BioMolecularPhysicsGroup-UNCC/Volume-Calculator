////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_FlatFileReader_hpp
#define bmpg_uncc_edu_util_FlatFileReader_hpp

#include <istream>
#include <string>
#include <cstdio>

namespace bmpg_uncc_edu {
namespace util {

/** FlatFileReader is a line oriented stream reader that erases
 *  leading and trailing whitespace and trailing comments.
*/
class FlatFileReader
{
public:
	FlatFileReader(std::istream & is, const std::string & comment_chars);
	~FlatFileReader();

	std::string get_line();
	operator bool() { return (is_ && is_.peek() != EOF); }
	size_t line_number() const { return line_number_; }

private:
	std::istream & is_;
	std::string comment_chars_;
	std::string whitespace_chars_;
	size_t line_number_;
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif


