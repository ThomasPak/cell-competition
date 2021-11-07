#ifndef DICTIONARY_HPP_
#define DICTIONARY_HPP_

#include <string>
#include <map>

/** Declaration of Dictionary alias.
 *
 * This file also contains a function to create a dictionary from an input
 * string.  The input string is expected to be formatted as newline-separated
 * strings of the form
 *
 * ```
 * KEY = VALUE
 * ```
 *
 * Whitespace surrounding the `=` symbol is ignored.
 */

/** Dictionary is a mapping from strings to strings, i.e.
 * std::map<std::string, std::string>.
 */
using Dictionary = std::map<std::string, std::string>;

/** Create dictionary from input string.
 *
 * The input string is expected to be formatted as newline-separated
 * strings of the form
 *
 * ```
 * KEY = VALUE
 * ```
 *
 * Whitespace surrounding the `=` symbol is ignored.
 *
 * @param input_string  input string.
 * @return created dictionary.
 */
Dictionary create_dictionary(const std::string& input_string);

bool is_whitespace(const char letter);
void remove_trailing_whitespace(std::string& str);
void remove_leading_whitespace(std::string& str);
void remove_leading_trailing_whitespace(std::string& str);

#endif // DICTIONARY_HPP_
