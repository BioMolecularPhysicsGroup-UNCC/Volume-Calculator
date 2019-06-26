////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_logger_hpp
#define bmpg_uncc_edu_util_logger_logger_hpp

#include <string>

namespace bmpg_uncc_edu {
namespace util {
namespace logger {
/**
 * Logger is a base class for logging messages, based on priority level.
 */
class Logger
{
	public:
		typedef enum { DEBUG, INFO, WARN, ERROR, CRITICAL } LOG_LEVEL;

	public:
		virtual std::ostream & stream();
		virtual void set_stream(std::ostream * os);
		virtual void set_log_level(LOG_LEVEL level);
		virtual void debug(const std::string & msg);
		virtual void info(const std::string & msg);
		virtual void warn(const std::string & msg);
		virtual void error(const std::string & msg);
		virtual void critical(const std::string & msg);
		virtual void log(const std::string & msg, LOG_LEVEL level) = 0;

	protected:
		Logger();
		virtual ~Logger();
		
	protected:
		LOG_LEVEL log_level_;
		std::ostream * os_;
};


} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

