/**
 * @file string_util.cpp
 *
 * Defines methods useful for formatting output.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#include "string_util.hpp"

using namespace mips;
using namespace mips::util;
using namespace std;

//! A utility function that replaces all all newlines with a number of spaces
//! depending on the indentation level.
string mips::util::Indent(string input)
{
  // Tab the first line.
  input.insert(0, 1, ' ');
  input.insert(0, 1, ' ');

  // Get the character sequence to replace all newline characters.
  std::string tabbedNewline("\n  ");

  // Replace all newline characters with the precomputed character sequence.
  size_t startPos = 0;
  while ((startPos = input.find("\n", startPos)) != string::npos)
  {
    input.replace(startPos, 1, tabbedNewline);
    startPos += tabbedNewline.length();
  }

  return input;
}
