////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_logger_cpp
#define bmpg_uncc_edu_util_logger_logger_cpp

#include <iostream>
#include <bmpg_uncc_edu/util/logger/Logger.hpp>

namespace bmpg_uncc_edu {
namespace util {
namespace logger {

Logger::Logger()
{
	log_level_ = INFO;
	os_ = &std::cout;
}

Logger::~Logger()
{
}

inline std::ostream & Logger::stream()
{
	return (*os_);
}

inline void Logger::set_stream(std::ostream * os)
{
	os_ = os;
}

inline void Logger::set_log_level(LOG_LEVEL level)
{
	log_level_ = level;
}

inline void Logger::debug(const std::string & msg)
{
	log(msg, DEBUG);
}

inline void Logger::info(const std::string & msg)
{
	log(msg, INFO);
}

inline void Logger::warn(const std::string & msg)
{
	log(msg, WARN);
}

inline void Logger::error(const std::string & msg)
{
	log(msg, ERROR);
}

inline void Logger::critical(const std::string & msg)
{
	log(msg, CRITICAL);
}


} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

