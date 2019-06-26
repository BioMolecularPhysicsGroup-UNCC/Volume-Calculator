////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_FlatFileReader_cpp
#define bmpg_uncc_edu_util_FlatFileReader_cpp

#include <bmpg_uncc_edu/util/FlatFileReader.hpp>

namespace bmpg_uncc_edu {
namespace util {

FlatFileReader::FlatFileReader(std::istream & is, const std::string & comment_chars) :
	is_(is), comment_chars_(comment_chars), whitespace_chars_(" \t\r\n\v\f"), line_number_(0)
{
}

FlatFileReader::~FlatFileReader()
{
}

std::string FlatFileReader::get_line()
{
	using namespace std;
	
	static const size_t buf_size = 1024;
	static char buffer[buf_size];
	size_t i;
	std::streamsize n;
	std::string line("");

	while (is_) {
		is_.getline(buffer, buf_size);
		++line_number_;
		n = is_.gcount();
		if (n <= 0)
			break;
		line = buffer;
		
		// Skip blank lines.  If not blank, erase leading whitespace.
		if ((i = line.find_first_not_of(whitespace_chars_)) == string::npos)
			continue;
		line.erase(0,i);
		// Erase trailing comments, if any
		if ((i = line.find_first_of(comment_chars_)) != string::npos)
			line.erase(i, line.length()-i);
		// Strip trailing whitespace
		if ((i = line.find_last_not_of(whitespace_chars_)) != string::npos)
			line.erase(i+1, line.length()-i-1);
		// Only return the line if it's not empty
		if (line.length() == 0)
			continue;
		else
			break;
	}
	
	return line;
}

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif

