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

#ifndef RG_IO_IOSTREAM_H
#define RG_IO_IOSTREAM_H

#include <util/config.h>

#include <util/io/manipulator.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <string>

namespace rg {

using std::string;

// ------------------------------------------------------------------------------------------------
// -- internal types and defaults for our string manipulators
// ------------------------------------------------------------------------------------------------	

namespace io { // hide 

/** Manipulators for stream output */ 
enum ManipulatorType { 
	openlist, /**< start of a list */
	closelist, /**< end of list */
	openitem, /**< start of an item */
	closeitem, /**< end of an item */
	delitem /**< delimiter between items */
};

/** total number of delimiters */
static const int NO_MANIPULATORS = delitem + 1;

static const string MANIPULATOR_DEFAULTS[NO_MANIPULATORS] = {
	"[", // openlist 
	"]", // closelist
	"",  // openitem
	"",  // closeitem
	"; " // delitem
};

} // namespace detail


// ------------------------------------------------------------------------------------------------
// -- internal functions
// ------------------------------------------------------------------------------------------------	

namespace io { // hide

/** Determines the index in the internal extensible array that corresponds to the given manipulator. 
 *
 * @param m type of manipulator
 * @return index of manipulator in the internal extensible array  
 */ 
inline int getHandle(ManipulatorType m) {
	static const int handles[NO_MANIPULATORS] = {
		allocManipulator(MANIPULATOR_DEFAULTS[openlist]),
		allocManipulator(MANIPULATOR_DEFAULTS[closelist]),
		allocManipulator(MANIPULATOR_DEFAULTS[openitem]),
		allocManipulator(MANIPULATOR_DEFAULTS[closeitem]),
		allocManipulator(MANIPULATOR_DEFAULTS[delitem])
	};
	
	return handles[m];
};

} // namespace io 

// ------------------------------------------------------------------------------------------------
// -- utility functions for setting manipulators
// ------------------------------------------------------------------------------------------------	

/** Set the string printed by StringUtils before the start of a list. */
inline Manipulator<string> set_openlist(const string& str) {
   return Manipulator<string>(io::getHandle(io::openlist), str);
};

/** Set the string printed by StringUtils before the end of a list. */
inline Manipulator<string> set_closelist(const string& str) {
   return Manipulator<string>(io::getHandle(io::closelist), str);
};

/** Set the string printed by StringUtils before the start of an item. */
inline Manipulator<string> set_openitem(const string& str) {
   return Manipulator<string>(io::getHandle(io::openitem), str);
};
	
/** Set the string printed by StringUtils before the end of an item. */
inline Manipulator<string> set_closeitem(const string& str) {
   return Manipulator<string>(io::getHandle(io::closeitem), str);
};

/** Set the delimiter printed by StringUtils in between items. */
inline Manipulator<string> set_delitem(const string& str) {
   return Manipulator<string>(io::getHandle(io::delitem), str);
};	


// ------------------------------------------------------------------------------------------------
// -- overload methods for stream output  
// ------------------------------------------------------------------------------------------------	

/** Output items of a sequence */
template<typename Iter, typename CharT, typename Traits>
void outputSequence(std::basic_ostream<CharT, Traits>& out, Iter begin, Iter end)
{
	using io::openitem;
	using io::delitem;
	using io::closeitem;
	using io::getHandle;	
	
	if (begin != end) {
		out << getManipulator<string>(out, getHandle(openitem)) 
			<< *begin 
			<< getManipulator<string>(out, getHandle(closeitem));
		++begin;
	}
	for(; begin != end; ++begin) {
		out << getManipulator<string>(out, getHandle(delitem))
			<< getManipulator<string>(out, getHandle(openitem)) 
			<< *begin 
			<< getManipulator<string>(out, getHandle(closeitem));
	};
}

/** Overload stream output for vectors */ 
template <typename T, typename CharT, typename Traits> 
std::basic_ostream<CharT, Traits>& operator<<(
			std::basic_ostream<CharT, Traits>& out, const std::vector<T>& v) 
{
	using io::openlist;
	using io::closelist;
	using io::getHandle;
	
	out << getManipulator<string>(out, getHandle(openlist));
	outputSequence(out, v.begin(), v.end());
	out << getManipulator<string>(out, getHandle(closelist));
	return out;
};

/** Overload stream output for pairs */
template <typename T1, typename T2, typename CharT, typename Traits>
std::basic_ostream<CharT, Traits>& operator<<(
			std::basic_ostream<CharT, Traits>& out, const std::pair<T1,T2>& v)
{
	using io::openlist;
	using io::openitem;
	using io::delitem;
	using io::closeitem;
	using io::closelist;
	using io::getHandle;

	out << getManipulator<string>(out, getHandle(openlist))
		<< getManipulator<string>(out, getHandle(openitem))
		<< v.first
		<< getManipulator<string>(out, getHandle(closeitem))
		<< getManipulator<string>(out, getHandle(delitem))
		<< getManipulator<string>(out, getHandle(openitem))
		<< v.second
		<< getManipulator<string>(out, getHandle(closeitem))
		<< getManipulator<string>(out, getHandle(closelist));

	return out;
};

} // namespace rg


#endif
