/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_Element_cpp
#define bmpg_uncc_edu_chemistry_Element_cpp

#include <iomanip>
#include <sstream>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/util/StringUtil.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {

size_t Element::Z() const
{
	return Z_;
}

std::ostream & operator<<(std::ostream & os, const Element & e)
{
	using namespace std;
	const char delim = ',';
	stringstream sstr;
	sstr << setw(3) << e.Z_ << delim;
	sstr << setw(3) << e.A_ << delim;
	sstr << setw(2) << e.symbol_ << delim;
	sstr << setw(15) << e.name_ << delim;
	sstr << scientific << setw(12) << setprecision(6) << e.chemical_weight_ << delim;
	sstr << scientific << setw(12) << setprecision(6) << e.atomic_mass_ << delim;
	sstr << scientific << setw(12) << setprecision(6) << e.abundance_ << delim;
        sstr << scientific << setw(12) << setprecision(6) << e.electronegativity_ << delim;
	sstr << setw(1) << e.radioactive_ << delim;
	if (e.radioactive_) {
		sstr << scientific << setw(12) << setprecision(6) << e.half_life_ << delim;
		for (size_t i = 0; i < e.decay_modes_.size(); i++) {
			Element::DECAY_MODE mode = e.decay_modes_[i];
			if (mode == Element::ALPHA_DECAY)
				sstr << "A:";
			else if (mode == Element::BETA_MINUS_DECAY)
				sstr << "B-:";
			else if (mode == Element::BETA_PLUS_DECAY)
				sstr << "B+:";
			else if (mode == Element::ELECTRON_CAPTURE_DECAY)
				sstr << "EC:";
		}
	}
	else
		sstr << delim << delim;  // indicate missing fields
	
	string s(sstr.str());
	char last_ch = s.at(s.length() - 1);
	if (!s.empty() && (last_ch == ':' || last_ch == delim)) // do not reverse order of comparison
		s=s.substr(0,s.length()-1); // eliminate trailing separator
	
	return (os << s); // Write the formatted string to the output stream and return
}

std::istream & operator>>(std::istream & is, Element & e)
{
	using namespace std;
	using namespace bmpg_uncc_edu::util;
	using namespace bmpg_uncc_edu::util::logger;
	Logger * logger = LoggerFactory::default_logger();

	std::string line;
	getline(is, line, '\n');
	std::vector<std::string> tokens = StringUtil::tokenize(line,true,",");
	if (tokens.size() != 10) {
		logger->warn("operator>>(istream &, Element &): (Z,A) pair doesn't have exactly 10 fields - skipping.");
		return is;
	}
	std::stringstream sbuf;
	for (size_t i = 0; i < tokens.size(); i++)
		sbuf << tokens[i] << '\n';
	std::stringstream sstr(sbuf.str());
	sstr >> e.Z_;//e._Z = StringUtil::str_to_size_t(tokens[0]);
	sstr >> e.A_;//e._A = StringUtil::str_to_size_t(tokens[1]);
	e.N_ = (e.A_ - e.Z_);
	sstr >> e.symbol_;//e._symbol = tokens[2];
	sstr >> e.name_; //e._name = tokens[3];
	sstr >> e.chemical_weight_; //e._chemical_weight = StringUtil::str_to_double(tokens[4]);
	sstr >> e.atomic_mass_;//e._atomic_mass = StringUtil::str_to_double(tokens[5]);
	sstr >> e.abundance_; //e._abundance = StringUtil::str_to_double(tokens[6]);
        sstr >> e.electronegativity_;
	sstr >> e.radioactive_; //e._radioactive = StringUtil::str_to_bool(tokens[7]);
	if (e.radioactive_) {
		sstr >> scientific >> e.half_life_; // e._half_life = StringUtil::str_to_double(tokens[8]);
		string mode_str;
		sstr >> mode_str;
		vector<string> mode_str_vec = StringUtil::tokenize(mode_str, false, ":");
		vector<Element::DECAY_MODE> mode_vec(mode_str_vec.size());
		for (size_t i = 0; i < mode_str_vec.size(); i++) {
			if (mode_str_vec[i] == "A")
				mode_vec.at(i) = Element::ALPHA_DECAY;
			else if (mode_str_vec[i] == "B-")
				mode_vec.at(i) = Element::BETA_MINUS_DECAY;
			else if (mode_str_vec[i] == "B+")
				mode_vec.at(i) = Element::BETA_PLUS_DECAY;
			else if (mode_str_vec[i] == "EC")
				mode_vec.at(i) = Element::ELECTRON_CAPTURE_DECAY;
			else
				logger->warn("operator>>(istream &, Element &): Unrecognized decay mode - skipping");
		}
		e.decay_modes_ = mode_vec;
	}
       	
	return is;
}

} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

