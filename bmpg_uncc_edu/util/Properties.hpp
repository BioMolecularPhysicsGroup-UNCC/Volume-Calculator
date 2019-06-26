////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_util_Properties_hpp
#define bmpg_uncc_edu_util_Properties_hpp

#include <string>
#include <map>

namespace bmpg_uncc_edu {
namespace util {

    /**
     * Properties splits a configuration text file into key/value pairs, removing white space 
     * and comments as appropriate.  The values are stored for retrieval by key.
     */
    
    
class Properties
{
	public:
		static const std::string NOT_FOUND;

	public:
		Properties();
		Properties(const Properties &);
		~Properties();

		// Loads the properties map from the specified input stream
		void load(std::istream & is, std::string comment_chars = "");
		// Writes the properties file to the specified output stream
		void list(std::ostream & os) const;
		// Returns true if the key exists, false otherwise
		bool has_key(std::string key) const;
		// Returns the value of the specified key or Properties::NOT_FOUND if not present
		std::string property(std::string key) const;
		// Returns the value of the specified key or a default value if not present
		std::string property(std::string key, const std::string & defaultValue) const;
		// Sets or inserts the key to have the specified value.
		void set_property(std::string key, std::string value);
		// Returns the underlying map used to store the key=value pairs
		const std::map<std::string, std::string> & map() const;
		
		std::string comment_characters() const { return comment_chars_; }
		void set_comment_characters(const std::string & cchars) { comment_chars_ = cchars; }

		friend std::istream & operator>>(std::istream & is, Properties & p);
		friend std::ostream & operator<<(std::ostream & os, const Properties & p);
		
	private:
		std::map<std::string, std::string> map_;
		std::string comment_chars_;
};

} // namespace bmpg_uncc_edu::util
} // namespace bmpg_uncc_edu

#endif


