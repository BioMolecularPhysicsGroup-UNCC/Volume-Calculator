////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_defaultlogger_cpp
#define bmpg_uncc_edu_util_logger_defaultlogger_cpp

#include <ctime>
#include <sstream>
#include <bmpg_uncc_edu/util/logger/DefaultLogger.hpp>

namespace bmpg_uncc_edu {
namespace util {
namespace logger {

DefaultLogger * DefaultLogger::instance_ = NULL;

DefaultLogger::DefaultLogger() : Logger()
{
	debug_count_ = info_count_ = warn_count_ = error_count_ = critical_count_ = 0;
}

DefaultLogger::~DefaultLogger()
{
}

DefaultLogger * DefaultLogger::instance()
{
	if (instance_ == NULL)
		instance_ = new DefaultLogger();
	return instance_;
}

void DefaultLogger::log(const std::string & msg, LOG_LEVEL level)
{
	// Do nothing if the log level is set to ignore this message
	if (level < log_level_)
		return;
	
	std::stringstream sstr;
	
	if (level == DEBUG) {
		sstr << "DEBUG : ";
		debug_count_++;
	}
	else if (level == INFO) {
		sstr << "INFO : ";
		info_count_++;
	}
	else if (level == WARN) {
		sstr << "WARNING : ";
		warn_count_++;
	}
	else if (level == ERROR) {
		sstr << "ERROR : ";
		error_count_++;
	}
	else if (level == CRITICAL) {
		sstr << "CRITICAL : ";
		critical_count_++;
	}
	else
		sstr << "UNKNOWN : ";
	
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	
	std::string datetime_str(asctime(timeinfo));
	if (datetime_str.at(datetime_str.size()-1) == '\n')
		datetime_str.at(datetime_str.size()-1) = ' ';
	sstr << datetime_str << ": " << msg;
	
	(*os_) << sstr.str() << "\n";
}

inline int DefaultLogger::debug_count()
{
	return debug_count_;
}

inline int DefaultLogger::warn_count()
{
	return warn_count_;
}

inline int DefaultLogger::info_count()
{
	return info_count_;
}

inline int DefaultLogger::error_count()
{
	return error_count_;
}

inline int DefaultLogger::critical_count()
{
	return critical_count_;
}


} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

