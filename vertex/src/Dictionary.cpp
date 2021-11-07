#include <stdexcept>
#include <sstream>

#include "Exception.hpp"

#include "Dictionary.hpp"

bool is_whitespace(const char letter)
{
    return (letter == ' ')
        || (letter == '\t')
        || (letter == '\n');
}

void remove_trailing_whitespace(std::string& str)
{
    auto pos = str.cend();

    // Iterate backwards until non-whitespace character is found
    while (is_whitespace(*--pos));

    // pos is off by one
    pos++;

    // Erase trailing whitespace and return string
    str.erase(pos, str.cend());
}

void remove_leading_whitespace(std::string& str)
{
    auto pos = str.cbegin();

    // Iterate backwards until non-whitespace character is found
    while (is_whitespace(*pos))
    {
        pos++;
    }

    // Erase trailing whitespace and return string
    str.erase(str.cbegin(), pos);
}

void remove_leading_trailing_whitespace(std::string& str)
{
    remove_leading_whitespace(str);
    remove_trailing_whitespace(str);
}

Dictionary create_dictionary(const std::string& input)
{
    // Initialise dictionary
    Dictionary dict;

    std::stringstream ss(input);
    std::string line;
    while (std::getline(ss, line))
    {
        // Split line by '='
        std::stringstream line_ss(line);

        std::string key;
        std::getline(line_ss, key, '=');

        // Skip line if no '=' is found
        if (line_ss.eof())
        {
            continue;
        }

        // Because I skip when no '=' is found, and I am using a stringstream,
        // it should not fail
        if (!line_ss)
        {
            NEVER_REACHED;
        }

        std::string value;
        std::getline(line_ss, value, '=');

        // Remove leading and trailing whitespace
        remove_leading_trailing_whitespace(key);
        remove_leading_trailing_whitespace(value);

        if ((key.size() == 0) || (value.size() == 0))
        {
            std::string error_msg;
            error_msg += "\"";
            error_msg += line;
            error_msg += "\" is missing key or value";
            EXCEPTION(error_msg);
        }

        dict[key] = value;
    }

    return dict;
}
