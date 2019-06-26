////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_LoggerFactory_cpp
#define bmpg_uncc_edu_util_logger_LoggerFactory_cpp

#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/util/logger/DefaultLogger.hpp>

namespace bmpg_uncc_edu {
namespace util {		
namespace logger {
	
Logger * LoggerFactory::default_logger_ = DefaultLogger::instance();

Logger * LoggerFactory::default_logger()
{
	return default_logger_;
}

void LoggerFactory::set_default_logger(Logger * logger_)
{
	default_logger_ = logger_;
}


} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

