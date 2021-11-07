#ifndef TESTDICTIONARY_HPP_
#define TESTDICTIONARY_HPP_

#include <cxxtest/TestSuite.h>

#include "Dictionary.hpp"

class TestDictionary : public CxxTest::TestSuite
{
    public:

        void TestParse1()
        {
            // Input string
            std::string input = "colour = green\nshape = triangle";

            // Create dictionary
            Dictionary dict = create_dictionary(input);

            // Assertions
            TS_ASSERT_EQUALS(dict.at("colour"), "green");
            TS_ASSERT_EQUALS(dict.at("shape"), "triangle");
        }

        void TestParse2()
        {
            // Input string
            std::string input = "  colour = \tgreen\n\nshape =   triangle  ";

            // Create dictionary
            Dictionary dict = create_dictionary(input);

            // Assertions
            TS_ASSERT_EQUALS(dict.at("colour"), "green");
            TS_ASSERT_EQUALS(dict.at("shape"), "triangle");
        }

        void TestParse3()
        {
            // Input string
            std::string input = "colour = gre en\n shape = tri\tan gle  ";

            // Create dictionary
            Dictionary dict = create_dictionary(input);

            // Assertions
            TS_ASSERT_EQUALS(dict.at("colour"), "gre en");
            TS_ASSERT_EQUALS(dict.at("shape"), "tri\tan gle");
        }

        void TestParse4()
        {
            // Input string
            std::string input = "colour = green = red\n shape = triangle";

            // Create dictionary
            Dictionary dict = create_dictionary(input);

            // Assertions
            TS_ASSERT_EQUALS(dict.at("colour"), "green");
            TS_ASSERT_EQUALS(dict.at("shape"), "triangle");
        }

        void TestError1()
        {
            // Faulty input string
            std::string input = " = green\n shape = triangle";

            // Assertions
            TS_ASSERT_THROWS_ANYTHING(Dictionary dict = create_dictionary(input););
        }

        void TestError2()
        {
            // Faulty input string
            std::string input = "colour =\n shape = triangle";

            // Assertions
            TS_ASSERT_THROWS_ANYTHING(Dictionary dict = create_dictionary(input););
        }

        void TestError3()
        {
            // Input string
            std::string input = "colour = green\n = triangle";

            // Assertions
            TS_ASSERT_THROWS_ANYTHING(Dictionary dict = create_dictionary(input););
        }

        void TestError4()
        {
            // Input string
            std::string input = "colour = green\nshape =  \t";

            // Assertions
            TS_ASSERT_THROWS_ANYTHING(Dictionary dict = create_dictionary(input););
        }
};

#endif // TESTDICTIONARY_HPP_
