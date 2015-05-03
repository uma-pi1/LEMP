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

#ifndef RG_IO_CONVERT_H
#define RG_IO_CONVERT_H

#include <iostream>
#include <sstream>
#include <string>

#include <util/exception.h>

namespace rg {

inline double toDouble(std::string const& s, bool failIfLeftoverChars = true) {
	std::istringstream i(s);
	double x;
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c))) {
		RG_THROW(InvalidArgumentException, string("Cannot convert to double: ") + s);
	}
	return x;
}

}

#endif
