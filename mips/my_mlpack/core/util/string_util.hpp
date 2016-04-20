/**
 * @file string_util.hpp
 *
 * Declares methods that are useful for writing formatting output.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_CORE_STRING_UTIL_HPP
#define __MY_MLPACK_CORE_STRING_UTIL_HPP

#include <string>

namespace mips {
namespace util {

//! A utility function that replaces all all newlines with a number of spaces
//! depending on the indentation level.
std::string Indent(std::string input);

}; // namespace util
}; // namespace mlpack

#endif
