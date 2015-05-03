
//    Copyright 2015 Rainer Gemulla
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#ifndef RG_EXCEPTION_CLASS_H
#define RG_EXCEPTION_CLASS_H

#include <util/config.h>

#include <string>
#include <exception>
#include <stdexcept>
#include <vector>

/** Shortcut for throwing exceptions with error message, filename and line number.
 * 
 * @param e class name of exception
 * @param m error message 
 */
#define RG_THROW(e, m) throw e(m, __FILE__, __LINE__)

namespace rg { 

using std::pair;
using std::string;
using std::vector;

// ------------------------------------------------------------------------------------------------
// -- Exception base class        
// ------------------------------------------------------------------------------------------------	

/** The base class for all exceptions in this library. Provides a formatted output that
 * includes the filename and line where the exception occured. On some systems, it also includes
 * a stack trace. */
class Exception : public std::exception {
	/** The error message supplied when throwing the exception */
	const string& m_errorMsg;
	
	/** The filename where the exception occured  */
	const string& m_file;
	
	/** The line at which the exception occured */
	const int m_line;
	
	/** number of stack calls to skip in formatted the stack trace */
	const int m_skip;
	
	/** The stack trace */
	vector< pair<string,string> > m_stackTrace;

	/** Maximum size of the error message */
	static const size_t FORMATTED_MSG_SIZE=4096;

	/** The fully formated error message including the stack trace */
	char m_formattedMsg[FORMATTED_MSG_SIZE];

	/** 
	 * Initializes this object. Derives the stack trace and formats the message output when
	 * what() is called. 
	 */ 
	void init();

protected:
	/** Constructor.
	 * 
	 * @param errorMsg the error message
	 * @param file the filename where the exception occured
	 * @param line the line at which the exception occured
	 * @param skip number of stack calls to skip in the formatted stack trace 
	 *             (increase by one for every level in the subclass hierarchy)
	 */ 
	explicit Exception(const string& errorMsg, const string& file, int line, int skip) throw() 
		: m_errorMsg(errorMsg), m_file(file), m_line(line), m_skip(skip) 
	{
		init();
	};
	
public:
	/** Constructor.
	 * 
	 * @param errorMsg the error message
	 * @param file the filename where the exception occured
	 * @param line the line at which the exception occured
	 */ 
	explicit Exception(const string& errorMsg = "", const string& file = "", int line = -1) throw() 
		: m_errorMsg(errorMsg), m_file(file), m_line(line), m_skip(0) 
	{
		init();
	};

	/** Destructor */
	~Exception() throw() { }
	

	/** Return a formatted error message.
	 * 
	 * @return the formatted error message
	 */
	const char* what() const throw() {
		return m_formattedMsg; 
	};
	
	/** Determine the stack trace (GNU compiler only, else does nothing). The stack trace is 
	 * appended to the supplied vector. Each element of the vector is a pair of two strings: 
	 * (1) the name of the binary and (2) the name and argument types of the method within 
	 * the binary. Note that the first element of the stack trace corresponds to the call of 
	 * this method. 
	 * 
	 * @param stackTrace vector in which to write the stack trace
	 * @param return true when successful
	 */
	static const bool getStackTrace(vector< pair<string,string> >& stackTrace);
};


// ------------------------------------------------------------------------------------------------
// -- exceptions used within this library
// ------------------------------------------------------------------------------------------------

/** Exception thrown when an operation is unsupported. Used in conjunction with "fat interfaces". */
class UnsupportedOperationException : public Exception {
public:
	explicit UnsupportedOperationException(const string& msg = "", 
			const string& file = "", int line = -1) throw() 
			: Exception(msg, file, line, 1) {};
};

/** Exception thrown when an argument is invalid. */
class InvalidArgumentException : public Exception {
public:
	explicit InvalidArgumentException(const string& msg = "", 
			const string& file = "", int line = -1) throw() 
			: Exception(msg, file, line, 1) {};
};

/** Exception thrown when an index is out of range. */
class OutOfRangeException : public Exception {
public:
	explicit OutOfRangeException(const string& msg = "", 
			const string& file = "", int line = -1) throw() 
			: Exception(msg, file, line, 1) {};
};

/** Exception thrown when something is not implemented but thought to be implemented. This exception 
 * is used to mark positions where an implementation is required but I just did not have the time
 * or need to do so. */
class NotImplementedException : public Exception {
public:
	explicit NotImplementedException(const string& msg = "", 
			const string& file = "", int line = -1) throw() 
			: Exception(msg, file, line, 1) {};
};


/** Exception thrown when the internal state is invalid. */
class IllegalStateException : public Exception {
public:
	explicit IllegalStateException(const string& msg = "", 
			const string& file = "", int line = -1) throw() 
			: Exception(msg, file, line, 1) {};
};

/** Exception thrown when I/O failed. */
class IOException : public Exception {
public:
	explicit IOException(const string& msg = "",
			const string& file = "", int line = -1) throw()
			: Exception(msg, file, line, 1) {};
};

} // namespace

#endif
