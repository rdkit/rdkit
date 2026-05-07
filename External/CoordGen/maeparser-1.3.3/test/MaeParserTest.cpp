#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/test/unit_test.hpp>

#include "Buffer.hpp"
#include "MaeBlock.hpp"
#include "MaeConstants.hpp"
#include "MaeParser.hpp"

using namespace schrodinger;
using namespace schrodinger::mae;

using std::shared_ptr;

using boost::trim_copy;

static std::string get_string(IndexedBlockBuffer& indexed_block_buffer,
                              size_t index)
{
    const char* ptr = nullptr;
    size_t len = 0;
    indexed_block_buffer.getData(index, &ptr, &len);
    return std::string(ptr, len);
}

BOOST_AUTO_TEST_SUITE(MaeParserSuite)

BOOST_AUTO_TEST_CASE(OuterBlockBeginning)
{
    {
        std::stringstream ss("{");
        Buffer b(ss);
        b.load();
        std::string name = outer_block_beginning(b);
        BOOST_REQUIRE_EQUAL(name, "");
    }
    {
        std::stringstream ss("f_m_ct {");
        Buffer b(ss);
        b.load();
        std::string name = outer_block_beginning(b);
        BOOST_REQUIRE_EQUAL(name, CT_BLOCK);
    }
    {
        std::stringstream ss("f_m_ct{");
        Buffer b(ss);
        b.load();
        std::string name = outer_block_beginning(b);
        BOOST_REQUIRE_EQUAL(name, CT_BLOCK);
    }
    {
        // TODO: Check that this is actually the desired behavior; I'm not
        //       sure that underscores are allowed in block names.
        std::stringstream ss("f_m_ct_block{");
        Buffer b(ss);
        b.load();
        std::string name = outer_block_beginning(b);
        BOOST_REQUIRE_EQUAL(name, "f_m_ct_block");
    }
}

BOOST_AUTO_TEST_CASE(OuterBlockBeginErrors)
{
    {
        std::stringstream ss("b_m_ct {");
        Buffer b(ss);
        b.load();
        try {
            std::string name = outer_block_beginning(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 1: "
                                "Bad format for outer block name; "
                                "must be (f|p)_<author>_<name>.");
        }
    }
    {
        std::stringstream ss("f_m {");
        Buffer b(ss);
        b.load();
        try {
            std::string name = outer_block_beginning(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 4: "
                                "Bad format for outer block name; "
                                "must be (f|p)_<author>_<name>.");
        }
    }
    {
        std::stringstream ss("full_m_ct {");
        Buffer b(ss);
        b.load();
        try {
            std::string name = outer_block_beginning(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 2: "
                                "Bad format for outer block name; "
                                "must be (f|p)_<author>_<name>.");
        }
    }
    {
        std::stringstream ss("f_m_ct   b_m_foo");
        Buffer b(ss);
        b.load();
        try {
            std::string name = outer_block_beginning(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 10: "
                                "Missing '{' for outer block.");
        }
    }
}

BOOST_AUTO_TEST_CASE(BlockBeginning)
{
    {
        auto ss = std::make_shared<std::stringstream>("m_something {");
        int indexed = -2;
        MaeParser mp(ss);
        std::string name = mp.blockBeginning(&indexed);
        BOOST_REQUIRE_EQUAL(name, "m_something");
        BOOST_REQUIRE_EQUAL(indexed, -1);
    }
    {
        auto ss = std::make_shared<std::stringstream>("mmmm_block{");
        int indexed = -2;
        MaeParser mp(ss);
        std::string name = mp.blockBeginning(&indexed);
        BOOST_REQUIRE_EQUAL(name, "mmmm_block");
        BOOST_REQUIRE_EQUAL(indexed, -1);
    }
    {
        auto ss = std::make_shared<std::stringstream>("m_whatev[23]{");
        int indexed = -2;
        MaeParser mp(ss);
        std::string name = mp.blockBeginning(&indexed);
        BOOST_REQUIRE_EQUAL(name, "m_whatev");
        BOOST_REQUIRE_EQUAL(indexed, 23);
    }
}

BOOST_AUTO_TEST_CASE(BlockBeginningErrors)
{
    {
        try {
            auto ss = std::make_shared<std::stringstream>("");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(
                trim_copy(std::string(e.what())),
                "Line 1, column 1: "
                "Bad format for block name; must be <author>_<name>.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("m_block[integer]");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 9: "
                                "Unexpected character.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("m_block[33  ");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 13: "
                                "Bad block index; missing ']'.");
        }
    }
    {
        try {
            auto ss =
                std::make_shared<std::stringstream>("m_block[33]  s_m_foo");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 14: "
                                "Missing '{' for block.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("'bad_block");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(
                trim_copy(std::string(e.what())),
                "Line 1, column 1: "
                "Bad format for block name; must be <author>_<name>.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("mmmm_ ");
            int indexed = 0;
            MaeParser mp(ss);
            mp.blockBeginning(&indexed);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(
                trim_copy(std::string(e.what())),
                "Line 1, column 6: "
                "Bad format for block name; must be <author>_<name>.");
        }
    }
}

BOOST_AUTO_TEST_CASE(BlockBody)
{
    {
        auto ss =
            std::make_shared<std::stringstream>(" b_m_foo b_m_bar ::: 1 0 }");
        MaeParser mp(ss);
        auto bl = mp.blockBody(CT_BLOCK);
        BOOST_REQUIRE(bl->getBoolProperty("b_m_foo"));
        BOOST_REQUIRE(!bl->getBoolProperty("b_m_bar"));
    }
    {
        auto ss = std::make_shared<std::stringstream>(
            " b_m_foo b_m_bar s_m_foo r_m_foo i_m_foo ::: "
            " 1       0       svalue  3.1415  22 }");
        MaeParser mp(ss);
        auto bl = mp.blockBody(CT_BLOCK);
        BOOST_REQUIRE(bl->getBoolProperty("b_m_foo"));
        BOOST_REQUIRE(!bl->getBoolProperty("b_m_bar"));
        BOOST_REQUIRE_EQUAL(bl->getStringProperty("s_m_foo"), "svalue");
        double tolerance = std::numeric_limits<double>::epsilon();
        BOOST_REQUIRE_CLOSE(bl->getRealProperty("r_m_foo"), 3.1415, tolerance);
        BOOST_REQUIRE_EQUAL(bl->getIntProperty("i_m_foo"), 22);
    }
}

BOOST_AUTO_TEST_CASE(BlockBodyErrors)
{
    {
        auto ss = std::make_shared<std::stringstream>(
            " b_m_foo\n s_m_foo\n r_m_foo\n i_m_foo\n :::\n"
            " 1\n svalue\n 3.1415\n 22\n ");
        MaeParser mp(ss);
        try {
            mp.blockBody(CT_BLOCK);
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 10, column 2: "
                                "Missing '}' for block.");
        }
    }
}

BOOST_AUTO_TEST_CASE(Property)
{
    {
        auto ss = std::make_shared<std::stringstream>("b_m_foo ");
        MaeParser mp(ss);
        auto p = mp.property();
        BOOST_REQUIRE(*p == "b_m_foo");
    }
    {
        auto ss = std::make_shared<std::stringstream>("r_m_bar ");
        MaeParser mp(ss);
        auto p = mp.property();
        BOOST_REQUIRE(*p == "r_m_bar");
    }
    {
        auto ss = std::make_shared<std::stringstream>("b_st_1_2_3_4_R_5 ");
        MaeParser mp(ss);
        auto p = mp.property();
        BOOST_REQUIRE(*p == "b_st_1_2_3_4_R_5");
    }
    {
        auto ss = std::make_shared<std::stringstream>("s_author_name ");
        MaeParser mp(ss);
        auto p = mp.property();
        BOOST_REQUIRE(*p == "s_author_name");
    }
}

BOOST_AUTO_TEST_CASE(PropertyValueSeparator)
{
    {
        auto ss = std::make_shared<std::stringstream>(":::");
        MaeParser mp(ss);
        auto p = mp.property();
        BOOST_REQUIRE(p == nullptr);
    }
}

BOOST_AUTO_TEST_CASE(PropertyList)
{
    {
        std::vector<std::string> properties;
        auto ss = std::make_shared<std::stringstream>("b_m_foo s_j_bar :::");
        // Buffer size of 14 was chosen to reproduce a bug during development;
        // a buffer boundary in the name part of the key.
        MaeParser mp(ss, 14u);
        while (true) {
            auto p = mp.property();
            mp.whitespace();
            if (p == nullptr) {
                break;
            }
            properties.push_back(*p);
        }
        BOOST_REQUIRE_EQUAL(properties.size(), 2u);
        BOOST_REQUIRE_EQUAL(properties[0], "b_m_foo");
        BOOST_REQUIRE_EQUAL(properties[1], "s_j_bar");
    }
    {
        std::vector<std::string> properties;
        auto ss = std::make_shared<std::stringstream>("b_m_foo s_j_ :::");
        MaeParser mp(ss);
        try {
            while (true) {
                auto p = mp.property();
                mp.whitespace();
                if (p == nullptr) {
                    break;
                }
                properties.push_back(*p);
            }
            BOOST_FAIL("Exception expected.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 13: "
                                "Bad format for property; "
                                "must be (b|i|r|s)_<author>_<name>.");
        }
    }
}

BOOST_AUTO_TEST_CASE(PropertyErrors)
{
    {
        try {
            auto ss = std::make_shared<std::stringstream>("bo_m_foo ");
            MaeParser mp(ss);
            mp.property();
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 2: "
                                "Bad format for property; "
                                "must be (b|i|r|s)_<author>_<name>.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("x_m_foo ");
            MaeParser mp(ss);
            mp.property();
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 1: "
                                "Bad format for property; "
                                "must be (b|i|r|s)_<author>_<name>.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("s_m_");
            MaeParser mp(ss);
            mp.property();
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 5: "
                                "Bad format for property; "
                                "must be (b|i|r|s)_<author>_<name>.");
        }
    }
    {
        try {
            auto ss = std::make_shared<std::stringstream>("s_m_ ");
            MaeParser mp(ss);
            mp.property();
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 5: "
                                "Bad format for property; "
                                "must be (b|i|r|s)_<author>_<name>.");
        }
    }
}

BOOST_AUTO_TEST_CASE(Integer)
{
    {
        std::stringstream ss("1234");
        Buffer b(ss);
        BOOST_REQUIRE_EQUAL(parse_value<int>(b), 1234);
    }
    {
        std::stringstream ss("-1234");
        Buffer b(ss);
        BOOST_REQUIRE_EQUAL(parse_value<int>(b), -1234);
    }
    {
        std::stringstream ss("2147483647");
        Buffer b(ss);
        BOOST_REQUIRE_EQUAL(parse_value<int>(b), 2147483647);
    }
    {
        // There is a bug in some VS editions that raises warning C4146
        // when assigning -2147483648 to an int in code..
        const int reference = std::numeric_limits<int>::min();
        std::stringstream ss;
        ss << reference;
        Buffer b(ss);
        BOOST_REQUIRE_EQUAL(parse_value<int>(b), reference);
    }
}

BOOST_AUTO_TEST_CASE(IntegerErrors)
{
    {
        std::stringstream ss("12-34");
        Buffer b(ss);
        try {
            parse_value<int>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 3: Unexpected '-'.");
        }
    }

    {
        std::stringstream ss("-12-34");
        Buffer b(ss);
        try {
            parse_value<int>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 4: Unexpected '-'.");
        }
    }

    {
        std::stringstream ss("\n\n]");
        Buffer b(ss);
        whitespace(b);
        try {
            parse_value<int>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 3, column 1: Missing integer.");
        }
    }

    {
        std::stringstream ss("\n\n123*]");
        Buffer b(ss);
        whitespace(b);
        try {
            parse_value<int>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 3, column 4: Unexpected character.");
        }
    }
}

BOOST_AUTO_TEST_CASE(Real)
{
    double tolerance = std::numeric_limits<double>::epsilon();

    {
        std::stringstream ss("-2.3 ");
        Buffer b(ss);
        BOOST_REQUIRE_CLOSE(parse_value<double>(b), -2.3, tolerance);
    }
    {
        std::stringstream ss("-24.3");
        Buffer b(ss);
        BOOST_REQUIRE_CLOSE(parse_value<double>(b), -24.3, tolerance);
    }
    {
        std::stringstream ss("-2.3e10 ");
        Buffer b(ss);
        BOOST_REQUIRE_CLOSE(parse_value<double>(b), -2.3e10, tolerance);
    }
    {
        std::stringstream ss("-2.3E10");
        Buffer b(ss);
        BOOST_REQUIRE_CLOSE(parse_value<double>(b), -2.3e10, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(RealErrors)
{
    {
        try {
            std::stringstream ss("");
            Buffer b(ss);
            parse_value<double>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 1: "
                                "Missing real.");
        }
    }
    {
        try {
            std::stringstream ss("-2.3{");
            Buffer b(ss);
            parse_value<double>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 5: "
                                "Unexpected character in real number.");
        }
    }
    {
        try {
            std::stringstream ss("\n -2.3. ");
            Buffer b(ss);
            whitespace(b);
            parse_value<double>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 2, column 6: "
                                "Bad real number.");
        }
    }
    {
        // Different versions of boost give different results due to
        // boost::spirit
        const std::string valid_error1("Line 2, column 2: Bad real number.");
        const std::string valid_error2("Line 2, column 4: Bad real number.");
        try {
            std::stringstream ss("\n -2EE3. ");
            Buffer b(ss);
            whitespace(b);
            parse_value<double>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_MESSAGE(
                trim_copy(std::string(e.what())) == valid_error1 ||
                    trim_copy(std::string(e.what())) == valid_error2,
                "Expected " << valid_error1 << " or " << valid_error2
                            << " but got " << e.what());
        }
    }
}

BOOST_AUTO_TEST_CASE(String)
{
    {
        std::stringstream ss("-2.3E10 ");
        Buffer b(ss);
        whitespace(b);
        std::string s = parse_value<std::string>(b);
        BOOST_REQUIRE_EQUAL(s, "-2.3E10");
    }
    {
        std::stringstream ss(R"(Q\ Z)");
        Buffer b(ss);
        whitespace(b);
        std::string s = parse_value<std::string>(b);
        BOOST_REQUIRE_EQUAL(s, "Q\\");
    }
    {
        std::stringstream ss(R"("Q\ Z")");
        Buffer b(ss);
        whitespace(b);
        std::string s = parse_value<std::string>(b);
        BOOST_REQUIRE_EQUAL(s, R"(Q Z)");
    }
    {
        std::stringstream ss(R"("a b c d e")");
        Buffer b(ss);
        whitespace(b);
        std::string s = parse_value<std::string>(b);
        BOOST_REQUIRE_EQUAL(s, "a b c d e");
    }
    {
        std::stringstream ss(" abcd"
                             "ef");
        Buffer b(ss, 5);
        whitespace(b);
        std::string s = parse_value<std::string>(b);
        BOOST_REQUIRE_EQUAL(s, "abcdef");
    }
}

BOOST_AUTO_TEST_CASE(StringErrors)
{
    {
        try {
            std::stringstream ss(R"("a b c d e)");
            Buffer b(ss);
            whitespace(b);
            parse_value<std::string>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 1, column 11: "
                                "Unterminated quoted string at EOF.");
        }
    }
}

BOOST_AUTO_TEST_CASE(Boolean)
{
    {
        std::stringstream ss(" 1");
        Buffer b(ss);
        whitespace(b);
        BOOST_REQUIRE_EQUAL(parse_value<BoolProperty>(b),
                            static_cast<BoolProperty>(true));
    }
    {
        std::stringstream ss("0 ");
        Buffer b(ss);
        whitespace(b);
        BOOST_REQUIRE_EQUAL(parse_value<BoolProperty>(b),
                            static_cast<BoolProperty>(false));
    }
}

BOOST_AUTO_TEST_CASE(BooleanErrors)
{
    {
        try {
            std::stringstream ss("\n\n\na");
            Buffer b(ss);
            whitespace(b);
            parse_value<BoolProperty>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 4, column 1: "
                                "Unexpected character for boolean value.");
        }
    }
    {
        try {
            std::stringstream ss("\t\n\n11");
            Buffer b(ss);
            whitespace(b);
            parse_value<BoolProperty>(b);
            BOOST_FAIL("Expected an exception.");
        } catch (read_exception& e) {
            BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                                "Line 3, column 2: "
                                "Unexpected character for boolean value.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(IndexedBlockBufferSuite)

BOOST_AUTO_TEST_CASE(TestIndexedValueCollectorInt)
{
    std::stringstream ss("   1 "
                         " <>  "
                         " 2");
    IndexedValueCollector<int> ivc("i_m_foo", 3);
    Buffer b(ss, 5);
    BOOST_REQUIRE(b.load());

    whitespace(b);
    ivc.parse(b);

    whitespace(b);
    ivc.parse(b);

    whitespace(b);
    ivc.parse(b);

    BOOST_REQUIRE_EQUAL(ivc.m_values[0], 1);
    BOOST_REQUIRE_EQUAL(ivc.m_values[1], 0);
    BOOST_REQUIRE_EQUAL(ivc.m_values[2], 2);
}

BOOST_AUTO_TEST_CASE(TestIndexedValueCollectorDouble)
{
    std::stringstream ss("   1 "
                         " <>  "
                         " 2.1");
    IndexedValueCollector<double> ivc("r_m_foo", 3);
    Buffer b(ss, 5);
    BOOST_REQUIRE(b.load());

    whitespace(b);
    ivc.parse(b);

    whitespace(b);
    ivc.parse(b);

    whitespace(b);
    ivc.parse(b);

    BOOST_REQUIRE_EQUAL(ivc.m_values[0], 1.0);
    BOOST_REQUIRE_EQUAL(ivc.m_values[1], 0.0);
    BOOST_REQUIRE_EQUAL(ivc.m_values[2], 2.1);
}

BOOST_AUTO_TEST_CASE(TestIndexedValueCollectorFunctors)
{
    std::stringstream ss("   1 "
                         " <>  "
                         " <>  "
                         " 2.1 "
                         " 3   "
                         " 0.1 ");
    Buffer b(ss, 5);
    IndexedValueCollector<double> real_ivc("r_m_foo", 3);
    IndexedValueCollector<int> int_ivc("i_m_foo", 3);

    std::vector<Parser*> parsers;
    parsers.push_back(&int_ivc);
    parsers.push_back(&real_ivc);

    BOOST_REQUIRE(b.load());
    for (int i = 0; i < 3; ++i) {
        for (auto parser : parsers) {
            whitespace(b);
            parser->parse(b);
        }
    }
    BOOST_CHECK_EQUAL(int_ivc.m_values[0], 1);
    BOOST_CHECK_EQUAL(int_ivc.m_values[1], 0);
    BOOST_CHECK_EQUAL(int_ivc.m_values[2], 3);
    BOOST_CHECK_EQUAL(real_ivc.m_values[0], 0.0);
    BOOST_CHECK_EQUAL(real_ivc.m_values[1], 2.1);
    BOOST_CHECK_EQUAL(real_ivc.m_values[2], 0.1);
}

// Split string.
BOOST_AUTO_TEST_CASE(Test0)
{
    std::stringstream ss("   1 "
                         " abc "
                         " ghij"
                         "k ");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghijk"));
}

// First string starts at column 1.
// Blank starting buffer.
// Blank buffer between strings.
// Blank buffer after all strings.
BOOST_AUTO_TEST_CASE(Test1)
{
    std::stringstream ss("   1 "
                         "     "
                         "abc  "
                         "     "
                         " ghij"
                         "k    "
                         "     ");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghijk"));
}

// Two strings that each fill a buffer and are separated by a blank buffer.
BOOST_AUTO_TEST_CASE(Test2)
{
    std::stringstream ss("   1 "
                         "abcde"
                         "     "
                         "fghij");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("f");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abcde"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("fghij"));
}

// Second string starts at buffer start, ends on next (and final) buffer end.
BOOST_AUTO_TEST_CASE(Test3)
{
    std::stringstream ss("   1 "
                         " abc "
                         "fghij"
                         "klmno");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("f");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("fghijklmno"));
}

// Second string starts at buffer start, ends on next (and final) buffer end.
BOOST_AUTO_TEST_CASE(Test3p1)
{
    std::stringstream ss("   1 "
                         " abc "
                         "fghij"
                         "klmno"
                         "pqrst");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("f");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("fghijklmnopqrst"));
}

// First string ends at buffer end.
// Second string spans an entire buffer in the middle.
BOOST_AUTO_TEST_CASE(Test4)
{
    std::stringstream ss("   1 "
                         "  abc"
                         " ghij"
                         "klmno"
                         "p  ");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghijklmnop"));
}

// Three separate strings with no split strings.
BOOST_AUTO_TEST_CASE(Test5)
{
    std::stringstream ss("   1 "
                         "     "
                         "  abc"
                         " ghi "
                         "jkl");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    ibb.addPropertyName("j");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghi"));
    z = get_string(ibb, 3);
    BOOST_REQUIRE_EQUAL(z, std::string("jkl"));
}

// Multi-buffer string that ends on a buffer end.
// Whitespace following last string doesn't end on buffer end.
BOOST_AUTO_TEST_CASE(Test6)
{
    std::stringstream ss("   1 "
                         "  a c"
                         "fghij"
                         " lmno"
                         "p  x ");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("c");
    ibb.addPropertyName("x");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("a"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("cfghij"));
    z = get_string(ibb, 3);
    BOOST_REQUIRE_EQUAL(z, std::string("lmnop"));
}

// Multiple strings in one buffer.
// Second buffer that's ignored.
BOOST_AUTO_TEST_CASE(Test7)
{
    std::stringstream ss("   1 "
                         "  abc ghi "
                         "xxxx   ");
    Buffer b(ss, 10);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("abc"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghi"));
}

// Unexpected EOF due to missing row number.
BOOST_AUTO_TEST_CASE(Test8)
{
    std::stringstream ss("  abc ghi ");
    Buffer b(ss, 10);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("a");
    ibb.addPropertyName("g");
    try {
        ibb.parse(b);
        BOOST_FAIL("Expected an exception.");
    } catch (read_exception& e) {
        BOOST_REQUIRE_EQUAL(trim_copy(std::string(e.what())),
                            "Line 1, column 11: "
                            "Unexpected EOF in indexed block values.");
    }
}

BOOST_AUTO_TEST_CASE(TestQuotedStrings1)
{
    std::stringstream ss(R"(   1      "bc  "     ghijk         )");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("b");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("\"bc  \""));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghijk"));
}

BOOST_AUTO_TEST_CASE(TestQuotedStrings2)
{
    std::stringstream ss(R"(   1      "bc \""   ghijk         )");
    Buffer b(ss, 5);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("b");
    ibb.addPropertyName("g");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("\"bc \\\"\""));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("ghijk"));
}

BOOST_AUTO_TEST_CASE(TestOneCharValues)
{
    std::stringstream ss("   1  1  2 3 40 ");
    Buffer b(ss, 128);
    IndexedBlockBuffer ibb("m_test", 1);
    ibb.addPropertyName("1");
    ibb.addPropertyName("2");
    ibb.addPropertyName("3");
    ibb.addPropertyName("40");
    ibb.parse(b);

    std::string z;
    z = get_string(ibb, 1);
    BOOST_REQUIRE_EQUAL(z, std::string("1"));
    z = get_string(ibb, 2);
    BOOST_REQUIRE_EQUAL(z, std::string("2"));
    z = get_string(ibb, 3);
    BOOST_REQUIRE_EQUAL(z, std::string("3"));
    z = get_string(ibb, 4);
    BOOST_REQUIRE_EQUAL(z, std::string("40"));
}

BOOST_AUTO_TEST_SUITE_END()
