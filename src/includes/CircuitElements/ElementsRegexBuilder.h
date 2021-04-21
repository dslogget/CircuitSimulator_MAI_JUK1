#ifndef _ELEMENTSREGEXBUILDER_HPP_INC_
#define _ELEMENTSREGEXBUILDER_HPP_INC_
#include <string>
#include <sstream>
#include <regex>

// indentifier is the sequence of letters representing the component
// simplifiedMatching:
// a string detailing how many, and the order of size_ts, and Values, and strings,
// and when options start startAnchor "^" endAnchor "$"
/// @brief a helper function to aid in the construction of regexes for parsing
///        netlist files
///
/// @param indentifier The designator for the component. e.g. "R" for resistor
/// @param simplifiedMatching The string that contains the data for what we want
///                           matched. For example:\n
/// n = int/size_t\n
/// w = word (works for floats etc)\n
/// ? = everything after is optional\n
/// c = char\n
/// s = space\n
/// @param startAnchor What prepends the regex, defaults to "^".
/// @param endAnchor What appends the regex, defaults to "$".
///
/// @return The complete Regex.
std::regex
generateRegex(std::string indentifier, std::string simplifiedMatching,
              bool startAnchor = true, bool endAnchor = true);

#endif
