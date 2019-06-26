////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_logger_defaultlogger_hpp
#define bmpg_uncc_edu_util_logger_defaultlogger_hpp

#include <bmpg_uncc_edu/util/logger/Logger.hpp>

namespace bmpg_uncc_edu {
namespace util {
namespace logger {    

/**
 * DefaultLogger defines the tracking and display of messages, based on priority level.
 */    

class DefaultLogger : public Logger
{
	public:
		static DefaultLogger * instance();
		virtual void log(const std::string & msg, LOG_LEVEL level);
		virtual int debug_count();
		virtual int warn_count();
		virtual int info_count();
		virtual int error_count();
		virtual int critical_count();

	protected:
		DefaultLogger();
		DefaultLogger(const DefaultLogger & rhs);
		DefaultLogger & operator=(const DefaultLogger & rhs);
		virtual ~DefaultLogger();
	
	private:
		static DefaultLogger * instance_;
		int debug_count_;
		int warn_count_;
		int info_count_;
		int error_count_;
		int critical_count_;
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu::util::logger
} // namespace bmpg_uncc_edu

#endif

