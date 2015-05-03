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


#ifndef RG_IO_UTILS_H
#define RG_IO_UTILS_H

#include <util/config.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

namespace rg {

using std::string;

/** Appends a string to another but ensures that no overflow occurs; the string to be
 * appended is truncated when necessary. This method does not use any stream 
 * operations and thus does not throw any exceptions.
 * 
 * @param to the target string
 * @param pos the position where to append (typically the index of the null value that
 *            terminates the target string). The position is modified to refer to the
 *            new end of the string.
 * @param from source string (null terminated)
 * @param maxLength maximum length of the target string, including the trailing zero
 */ 
inline void safeAppend(char* to, int& pos, const char* from, int maxLength) throw() {
	int len = std::min((int)strlen(from), maxLength-pos-1);
	if (len > 0) {
		strncpy(&to[pos], from, len);
		pos += len;
		to[pos] = 0;
	}
}

inline string underline(const string& s) {
	return s + "\n" + string(s.size(), '-') + "\n";	
}

template<typename T1>
inline std::string paste(const T1& a1) {
	std::stringstream ss;
	ss << a1;
	return ss.str();
}

template<typename T1, typename T2>
inline std::string paste(const T1& a1, const T2& a2) {
	std::stringstream ss;
	ss << a1 << a2;
	return ss.str();
}

template<typename T1, typename T2, typename T3>
inline std::string paste(const T1& a1, const T2& a2, const T3& a3) {
	std::stringstream ss;
	ss << a1 << a2 << a3;
	return ss.str();
}

template<typename T1, typename T2, typename T3, typename T4>
inline std::string paste(const T1& a1, const T2& a2, const T3& a3, const T4& a4) {
	std::stringstream ss;
	ss << a1 << a2 << a3 << a4;
	return ss.str();
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
inline std::string paste(const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5) {
	std::stringstream ss;
	ss << a1 << a2 << a3 << a4 << a5;
	return ss.str();
}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
inline std::string paste(const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6) {
	std::stringstream ss;
	ss << a1 << a2 << a3 << a4 << a5 << a6;
	return ss.str();
}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
inline std::string paste(const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T6& a7) {
	std::stringstream ss;
	ss << a1 << a2 << a3 << a4 << a5 << a6 << a7;
	return ss.str();
}


} // namespace

#endif
