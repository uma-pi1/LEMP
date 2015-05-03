
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

#include <util/exception/Exception.h>
#include <util/io.h>

#include <cstring>

// include GNU specific headers (for stacktrace) 
#ifdef __GNUC__
#include <execinfo.h>
#include <stdlib.h>
#include <cxxabi.h>
#include <stdio.h>
#endif


namespace rg {

// ------------------------------------------------------------------------------------------------
// -- Exception members  
// ------------------------------------------------------------------------------------------------	

void Exception::init() {
	// current position and length of string to be appended
	int pos = 0;
	
	// append message
	safeAppend(m_formattedMsg, pos, m_errorMsg.c_str(), FORMATTED_MSG_SIZE);
	safeAppend(m_formattedMsg, pos, "\n", FORMATTED_MSG_SIZE);

	// append file and line, if present
	if (m_line >= 0) {
		safeAppend(m_formattedMsg, pos, "    in():  ", FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, m_file.c_str(), FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, "(", FORMATTED_MSG_SIZE);
		
		// append linenumber
		char sline[20];
		sprintf(sline, "%d", m_line);
		safeAppend(m_formattedMsg, pos, sline, FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, ")\n", FORMATTED_MSG_SIZE);
	}
	
	// determine stack trace
	getStackTrace(m_stackTrace);
	
	// append stack trace (but skip the methods of this class)
	int i0 = 2 + m_skip;
	for (unsigned int i=i0; i<m_stackTrace.size(); i++) {
		// append spaces
		if (i == i0) {
			safeAppend(m_formattedMsg, pos, "    at():  ", FORMATTED_MSG_SIZE);
		} else {
			safeAppend(m_formattedMsg, pos, "           ", FORMATTED_MSG_SIZE);
		}
		
		safeAppend(m_formattedMsg, pos, m_stackTrace[i].first.c_str(), FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, ":", FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, m_stackTrace[i].second.c_str(), FORMATTED_MSG_SIZE);
		safeAppend(m_formattedMsg, pos, "\n", FORMATTED_MSG_SIZE);
	}
};

const bool Exception::getStackTrace(vector< pair<string,string> >& stackTrace) {
#ifdef __GNUC__
	const size_t max_depth = 100;
    size_t stack_depth;
    void *stack_addrs[max_depth];
    char **stack_strings;

    stack_depth = backtrace(stack_addrs, max_depth);
    stack_strings = backtrace_symbols(stack_addrs, stack_depth);

	for (size_t i = 1; i < stack_depth; i++) {
		size_t sz = 200; // just a guess, template names will go much wider
	    char *function = static_cast<char *>(malloc(sz));
	    char *begin = 0, *end = 0;
	    // find the parentheses and address offset surrounding the mangled name
	    for (char *j = stack_strings[i]; *j; ++j) {
	        if (*j == '(') {
	            begin = j;
	        }
	        else if (*j == '+') {
	            end = j;
	        }
	    }
	    if (begin && end) {
	        *begin++ = 0;
	        *end = 0;
	        // found our mangled name, now in [begin, end)
	
	        int status;
	        char *ret = abi::__cxa_demangle(begin, function, &sz, &status);
	        if (ret) {
	            // return value may be a realloc() of the input
	            function = ret;
	        } else {
	            // demangling failed, just pretend it's a C function with no args
	            std::strncpy(function, begin, sz);
	            std::strncat(function, "()", sz);
	            function[sz-1] = 0;
	        }
	        pair<string,string> p(stack_strings[i], function);
	        stackTrace.push_back(p);		        
	    }
	    else
	    {
	    	pair<string,string> p(stack_strings[i], "");
	        stackTrace.push_back(p);
	    }
	    free(function);
	}    
    free(stack_strings); // malloc()ed by backtrace_symbols
    return true;
#else
    return false; // not GNU
#endif
    
};

} // namespace
