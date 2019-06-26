/**
  Copyright (c) 2008 by Mike Fairchild
  @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
*/

#ifndef bmpg_uncc_edu_chemistry_library_PeriodicTable_cpp
#define bmpg_uncc_edu_chemistry_library_PeriodicTable_cpp

#include <memory>
#include <iostream>
#include <string>
#include <cctype>
#include <sstream>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>
#include <bmpg_uncc_edu/util/StringUtil.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
using namespace std;
using namespace::bmpg_uncc_edu::util::logger;

// static member initialization
PeriodicTable * PeriodicTable::instance_ = NULL;

// Private constructor (copy constructor, operator=, destructor are defaulted)
PeriodicTable::PeriodicTable()
{
}

PeriodicTable * PeriodicTable::instance()
{
	if (instance_ == NULL)
		instance_ = new PeriodicTable();
	return instance_;
}

const Element * const PeriodicTable::element(size_t Z,
					     size_t A) const
{
	using namespace bmpg_uncc_edu::util::logger;
	Logger * logger = LoggerFactory::default_logger();

	Element * e = NULL;
	double abundance_max = 0.0;
	table_t::const_iterator pos;
	for (pos = table_.begin(); pos != table_.end(); ++pos) {
		std::pair<size_t, size_t> za = pos->first;
		if (za.first == Z) {
			if (A != 0 && za.second == A)
				return pos->second;
			else if (A == 0) {
				if (pos->second->abundance() > abundance_max) {
					abundance_max = pos->second->abundance();
					e = pos->second;
				}
			}
		}
	}
	if (e == NULL){
		stringstream msg;
		msg << "PeriodicTable::element(size_t Z, size_t A): Attempt to find element failed; Z = " << Z << ", A = " << A << ".";
		logger->warn(msg.str());
	}
	return e;
}

const Element * const PeriodicTable::element(const std::string & symbol,
					     size_t A) const
{
	using namespace bmpg_uncc_edu::util;
	using namespace bmpg_uncc_edu::util::logger;
	Logger * logger = LoggerFactory::default_logger();
	const Element * e = NULL;

	// Make sure first letter is upper case, second letter is lower case, and
	// then chop off the rest.
	std::string sym(symbol);
	if (sym.length() >= 1)
		sym.at(0) = toupper(sym.at(0));
	if (sym.length() >= 2) {
		sym.at(1) = tolower(sym.at(1));
		sym = sym.substr(0,2); // we only want to search based on the first two characters
	}
	// Find the Z for the given symbol and then return element(Z,A)
	table_t::const_iterator pos;
	for (pos = table_.begin(); pos != table_.end(); ++pos) {
		if (pos->second->symbol() == sym) // case sensitive comparison
			e = element(pos->second->Z(), A);
	}
	if (e == NULL){
		stringstream msg;
		msg << "PeriodicTable::element(string symbol, size_t A): Attempt to find element failed; symbol = " << symbol << ", A = " << A << ".";
		logger->warn(msg.str());
	}
	return e;
}

std::istream & PeriodicTable::load(std::istream & is)
{
	return (is >> (*this));
}

std::istream & operator>>(std::istream & is, PeriodicTable & pt)
{
	using namespace std;
	using bmpg_uncc_edu::util::StringUtil;
	
	Logger * logger = LoggerFactory::default_logger();
	std::string line;
	while (is && !is.eof()) {
		getline(is, line, '\n');
		StringUtil::chop_whitespace(line);
		if (line.empty() || line.at(0) == '#') // do NOT reverse order of comparison
			continue; // Skip line if empty or comment
		std::stringstream sstr(line);

		auto_ptr<Element> elem(new Element());
		try {
			sstr >> *elem;
		}
		catch (exception & e) {
			logger->error(e.what());
		}
		size_t Z = elem->Z(), A = elem->A();
		if (!pt.table_.insert(std::make_pair(std::make_pair(Z,A),elem.get())).second)
			logger->warn("Duplicate (Z,A) entry found - skipping.\n");
		elem.release();
	}
	return is;
}

std::ostream & operator<<(std::ostream & os, const PeriodicTable & pt)
{
	PeriodicTable::table_t::const_iterator pos;
	for (pos = pt.table_.begin(); pos != pt.table_.end(); ++pos)
		os << *(pos->second) << '\n';
	return os;
}

}	//namespace bmpg_uncc_edu::chemistry::library
} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

