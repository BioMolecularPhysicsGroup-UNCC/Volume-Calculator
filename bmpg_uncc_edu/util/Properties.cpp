////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_Properties_cpp
#define bmpg_uncc_edu_util_Properties_cpp

#include <iostream>
#include <cctype>
#include <bmpg_uncc_edu/util/Properties.hpp>
#include <bmpg_uncc_edu/util/StringUtil.hpp>
#include <bmpg_uncc_edu/util/FlatFileReader.hpp>

namespace bmpg_uncc_edu {
namespace util {
using namespace std;

const std::string Properties::NOT_FOUND = "bmpg_uncc_edu::properties::NOT_FOUND";

Properties::Properties() : comment_chars_("#")
{
}

Properties::Properties(const Properties & copy) : comment_chars_("#")
{
	map_ = copy.map_;
}

Properties::~Properties()
{
	map_.clear();
}

void Properties::load(std::istream & is, std::string comment_chars)
{
	if (comment_chars == "")
		comment_chars = comment_chars_;

	std::string line;
	std::string key;
	std::string value;
	FlatFileReader reader(is, comment_chars_);
	while (reader) {
		line = reader.get_line();
                cout << line << endl;                                // DJJacobs  (to look)
		size_t pos = line.find_first_of('=');
		if (pos != std::string::npos) {
			key = line.substr(0,pos);
			value = line.substr(pos+1,line.length()-pos-1);
			
		}
		else {
			key = line;
			value = "";
		}
		
		// Remove leading and trailing whitespaces from key and value
		StringUtil::chop_whitespace(key);
		StringUtil::chop_whitespace(value);
		
                cout << "putative:  key = " << key << "  value = " << value << endl;
		// Skip this line if it's empty
		if (key.length() == 0)
			continue;

                cout << "actual:    key = " << key << "  value = " << value << endl;
		
		map_[key]=value;
	}
}

void Properties::list(std::ostream & os) const
{
	std::map<std::string, std::string>::const_iterator iter;
	for (iter = map_.begin(); iter != map_.end(); iter++)
		os << iter->first << "=" << iter->second << std::endl;
}

bool Properties::has_key(std::string key) const
{
	StringUtil::chop_whitespace(key);
	return (map_.find(key) != map_.end());
}

std::string Properties::property(std::string key) const
{
	return property(key, NOT_FOUND);
}

std::string Properties::property(std::string key, const std::string & defaultValue) const
{
	StringUtil::chop_whitespace(key);
	std::map<std::string,std::string>::const_iterator it;
	it = map_.find(key);
	if (it == map_.end())
		return defaultValue;
	else
		return it->second;
}

void Properties::set_property(std::string key, std::string value)
{
	StringUtil::chop_whitespace(key);
	StringUtil::chop_whitespace(value);
	map_[key]=value;
}

const std::map<std::string, std::string> & Properties::map() const
{
	return map_;
}

std::istream & operator>>(std::istream & is, Properties & p)
{
	p.load(is);
	return is;
}

std::ostream & operator<<(std::ostream & os, const Properties & p)
{
	p.list(os);
	return os;
}

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif


