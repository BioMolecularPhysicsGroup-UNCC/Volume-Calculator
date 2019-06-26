////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_StringUtil_hpp
#define bmpg_uncc_edu_util_StringUtil_hpp

#include <string>
#include <cstddef>
#include <vector>

namespace bmpg_uncc_edu {
namespace util {
    
    /**
     * StringUtil provides a variety of string formatting functions such as eliminating spaces, converting to other types
     * and tokenizing.
     */

class StringUtil
{
	public:
		const static std::string whitespaces_;
		
	public:
		static std::string & chop_whitespace(std::string &);
		static std::string & chop_trailing_whitespace(std::string &);
		static std::string & chop_leading_whitespace(std::string &);
		static std::vector<std::string> tokenize(const std::string & s, bool empty_tokens = false, const std::string delims = whitespaces_);
		static std::string pretty_format_bytes(double bytes);
		static bool str_to_bool(const std::string & str);
		static double str_to_double(const std::string & str);
		static size_t str_to_size_t(const std::string & str);
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif


