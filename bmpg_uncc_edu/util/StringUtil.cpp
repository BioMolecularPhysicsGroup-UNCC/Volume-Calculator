////////////////////////////////////////////////////////////////////////////////  
//  
// Copyright (c) 2008 by Mike Fairchild  
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com  
//  
////////////////////////////////////////////////////////////////////////////////  
  
#ifndef bmpg_uncc_edu_util_StringUtil_cpp  
#define bmpg_uncc_edu_util_StringUtil_cpp  
  
#include <cmath>  
#include <sstream>  
#include <limits>  
#include <cstdlib>  
#include <bmpg_uncc_edu/util/StringUtil.hpp>  
#include <bmpg_uncc_edu/util/Exception.hpp>  
  
namespace bmpg_uncc_edu {  
namespace util {  
  
const std::string StringUtil::whitespaces_ = " \t\f\v\n\r";  
  
std::string & StringUtil::chop_whitespace(std::string & str)  
{  
	chop_leading_whitespace(str);  
	chop_trailing_whitespace(str);  
	return str;  
}  
  
std::string & StringUtil::chop_leading_whitespace(std::string & str)  
{  
	size_t found;  
	found = str.find_first_not_of(whitespaces_);  
	if (found != std::string::npos)  
		str.erase(0,found); // get rid of leading whitespaces  
	else  
		str.clear(); // str is all whitespace  
	return str;  
}  
  
std::string & StringUtil::chop_trailing_whitespace(std::string & str)  
{  
	size_t found;  
	found = str.find_last_not_of(whitespaces_);  
	if (found != std::string::npos)  
		str.erase(found+1); // get rid of trailing whitespaces  
	else  
		str.clear();  // str is all whitespace  
	return str;  
}  
  
std::string StringUtil::pretty_format_bytes(double bytes)  
{  
	static const double k = 1024.0;  
	static const double mb = k*k;  
	static const double gb = k*mb;  
	bytes = fabs(bytes);  
	std::stringstream sstr;  
	if (bytes < k)  
		sstr << bytes << " bytes";  
	else if (bytes < mb)  
		sstr << (bytes/k) << " K";  
	else if (bytes < gb)  
		sstr << (bytes/mb) << " MB";  
	else  
		sstr << (bytes/gb) << " GB";  
	return sstr.str();  
}  
  
std::vector<std::string> StringUtil::tokenize(const std::string & s, bool empty_tokens, const std::string delims)  
{  
	std::vector<std::string> tokens;  
	size_t i = 0, j = 0;  
	if (s.empty()) // Short circuit if s is empty, else indexing below will fail (also more efficient)  
		return tokens;  
	if (delims.find_first_of(s.at(0)) != std::string::npos && empty_tokens)  
		tokens.push_back(""), ++i;  
	while ((j = s.find_first_of(delims, i)) != std::string::npos) {  
		if ((j - i) == 0 && empty_tokens)  
			tokens.push_back("");  
		else if ((j - i) >= 1)  
			tokens.push_back(s.substr(i,j-i));  
		i = j+1;  
	}  
	if ((j = delims.find_first_of(s.at(s.length() - 1))) != std::string::npos && empty_tokens)  
		tokens.push_back("");  
	else if (j == std::string::npos)  
		tokens.push_back(s.substr(i));  
  
	return tokens;  
}  
  
bool StringUtil::str_to_bool(const std::string & str)  
{  
	return (str_to_double(str) != 0.0);  
}  
  
double StringUtil::str_to_double(const std::string & str)  
{  
	return strtod(str.c_str(), NULL);  
}  
  
size_t StringUtil::str_to_size_t(const std::string & str)  
{  
	const char * s = str.c_str();   
	size_t result = std::numeric_limits<size_t>::max();  
	double d = strtod(s,NULL);  
	if (d < 0)  
		throw Exception("StringUtil::str_to_size_t(): String is negative.");  
	else if (round(d) != d)  
		throw Exception("StringUtil::str_to_size_t(): String is not an integer type.");  
	// FIXME - can a double exactly represent every number that a size_t can?  
	// I think so, but more thought should go into the next cast.  
	// There won't be problems with magnitudes, but possibly with the mantissa (precision)  
	result = static_cast<size_t>(d);  
	return result;  
}  
  
} // namespace bmpg_uncc_edu::util  
} // namespace bmpg_uncc_edu  
  
#endif  
  

