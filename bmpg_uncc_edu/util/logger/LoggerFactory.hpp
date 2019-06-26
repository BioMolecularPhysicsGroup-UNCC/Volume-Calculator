////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_LoggerFactory_hpp
#define bmpg_uncc_edu_util_logger_LoggerFactory_hpp

#include <bmpg_uncc_edu/util/logger/Logger.hpp>

namespace bmpg_uncc_edu {
namespace util {
namespace logger {
 
/**
 * LoggerFactory sets the default Logger class, used to log messages according to their priority 
 */
class LoggerFactory
{
	public:
		static Logger * default_logger();
		static void set_default_logger(Logger * logger_);

	private:
		// hide constructors, operator=, and destructor
		LoggerFactory();
		LoggerFactory(const LoggerFactory & rhs);
		LoggerFactory & operator=(const LoggerFactory & rhs);
		~LoggerFactory();
		
		static Logger * default_logger_;
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

