#include "CircuitElements/ElementsRegexBuilder.h"
// indentifier is the sequence of letters representing the component
// simplifiedMatching:
// a string detailing how many, and the order of size_ts, and Values, and strings, and when options start
// startAnchor "^"
// endAnchor "$"
//
// n = int/size_t
// w = word (works for floats etc)
// ? = everything after is optional
// c = char
// s = space
//
// Capacitor example:
//    ( "C", "n n w" )
// Resistor example:
//    ( "R", "n n w ? w" )
//
//

std::regex
generateRegex( std::string indentifier, std::string simplifiedMatching,
               bool startAnchor,
               bool endAnchor ) {
   constexpr char spaceRegex[] = R"(\s)";
   constexpr char startAnchorRegex[] = R"(^)";
   constexpr char endAnchorRegex[] = R"(\s?$)";
   constexpr char optionalSpaceRegex[] = R"(\s?)";

   constexpr char charRegex[] = R"((.))";
   constexpr char size_tRegex[] = R"((\d+?))";
   constexpr char emptyWordRegex[] = R"((.*?))";
   constexpr char wordRegex[] = R"((.+?))";
   constexpr char optionalCharRegex[] = R"((?:\s(.))?)";
   constexpr char optionalWordRegex[] = R"((?:\s(.+?))?)";
   constexpr char optionalSize_tRegex[] = R"((?:\s(\d+?))?)";

   std::string built( "" );

   if ( startAnchor ) {
      built += startAnchorRegex;
   }

   built += indentifier + emptyWordRegex;

   bool option = false;
   for ( char letter : simplifiedMatching ) {
      switch( letter ) {
         case 'n':
            if ( option ) {
               built += optionalSize_tRegex;
            } else {
               built += spaceRegex;
               built += size_tRegex;
            }
            break;
         case 'w':
            if ( option ) {
               built += optionalWordRegex;
            } else {
               built += spaceRegex;
               built += wordRegex;
            }
            break;
         case 'c':
            if ( option ) {
               built += optionalCharRegex;
            } else {
               built += spaceRegex;
               built += charRegex;
            }
            break;
         case 's':
            if ( option ) {
               built += optionalSpaceRegex;
            } else {
               built += spaceRegex;
            }
            break;
         case '?':
            option = true;
            break;;
         case ' ':
            break;
         default:
            break;
      }
   }

   if ( endAnchor ) {
      built += endAnchorRegex;
   }

   return std::regex( built.c_str() );
}
