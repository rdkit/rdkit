#include <cstdlib>
#include <cstring>

#include <map>
#include <sstream>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include "Buffer.hpp"

using namespace schrodinger;

BOOST_AUTO_TEST_SUITE(BufferSuite)

BOOST_AUTO_TEST_CASE(NewlineColumn0)
{
    std::string s("123456");
    std::stringstream ss(s);
    Buffer b(ss, s.size());
    BOOST_REQUIRE(b.load());

    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 0), 1u);
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 2), 3u);
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 5), 6u);
}

BOOST_AUTO_TEST_CASE(NewlineColumn1)
{
    std::string s("\n123456");
    std::stringstream ss(s);
    Buffer b(ss, s.size());
    BOOST_REQUIRE(b.load());

    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 0), 1u);
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 1), 1u);
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 3), 3u);
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin + 6), 6u);
}

BOOST_AUTO_TEST_CASE(NewlineColumnEmptyBuffer)
{
    std::string s;
    std::stringstream ss(s);
    Buffer b(ss, s.size());
    BOOST_REQUIRE(!b.load());
    BOOST_REQUIRE_EQUAL(b.getColumn(b.begin), 1u);
}

BOOST_AUTO_TEST_CASE(Loading)
{
    std::stringstream ss("123456");
    Buffer b(ss, 6);
    int i = 0;
    while (b.current < b.end || b.load()) {
        ++b.current;
        ++i;
    }
    BOOST_REQUIRE_EQUAL(i, 6);
}

BOOST_AUTO_TEST_CASE(GetBuffer)
{
    std::stringstream ss("123456123456");
    int buffer_len = 6;
    Buffer b(ss, buffer_len);
    for (int i = 0; i < buffer_len; i++) {
        if (!b.load()) {
            BOOST_FAIL("Unexpected end of \"file\".");
        }
        ++b.current;
    }
    if (!b.load()) {
        BOOST_FAIL("Unexpected end of \"file\".");
    }
    BOOST_REQUIRE_EQUAL(*b.current, '1');
}

BOOST_AUTO_TEST_CASE(Newline0)
{
    // Test newline as the first character of a new buffer load.
    std::stringstream ss("123456\n");
    int buffer_len = 6;
    Buffer b(ss, buffer_len);
    for (int i = 0; i < buffer_len; i++) {
        if (!b.load()) {
            BOOST_FAIL("Unexpected end of \"file\".");
        }
        ++b.current;
    }
    BOOST_REQUIRE_EQUAL(b.getColumn(), 7u);
    if (!b.load()) {
        BOOST_FAIL("Unexpected end of \"file\".");
    }
    BOOST_REQUIRE_EQUAL(*b.current, '\n');
    ++b.current;
    BOOST_REQUIRE_EQUAL(b.getColumn(), 1u);
}

BOOST_AUTO_TEST_CASE(GetPosition0)
{
    // Test expected column numbers for repeated buffer loads.
    std::stringstream ss("0123456\n12345x7\n\n123y5\n");
    //                    0......1.......2........3.
    //                    12345678 12345678 1 123456
    const unsigned int buffer_len = 7;
    const unsigned int stream_len = 23;
    Buffer b(ss, buffer_len);

    const int required_reads = 4;
    unsigned int column[stream_len] = {1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4,
                                       5, 6, 7, 8, 1, 1, 2, 3, 4, 5, 6};

    int stream_ix = 0;
    int load_count = 0;
    while (b.current < b.end || (b.load() && ++load_count)) {
        BOOST_REQUIRE_EQUAL(b.getColumn(), column[stream_ix]);
        ++b.current;
        ++stream_ix;
    }
    BOOST_REQUIRE_EQUAL(load_count, required_reads);
}

BOOST_AUTO_TEST_CASE(GetPosition1)
{
    // Test expected column numbers for repeated buffer loads with
    // the smallest possible buffer size.
    std::stringstream ss("0123456");

    int buffer_len = 1;
    Buffer b(ss, buffer_len);

    const int required_reads = 7;
    unsigned int column[required_reads] = {1, 2, 3, 4, 5, 6, 7};

    int load_count = 0;
    while (b.current < b.end || (b.load() && ++load_count)) {
        BOOST_REQUIRE_EQUAL(b.getColumn(), column[load_count - 1]);
        ++b.current;
    }
    BOOST_REQUIRE_EQUAL(load_count, required_reads);
}

BOOST_AUTO_TEST_CASE(NewlineCase0)
{
    // Test the corner case of a single newline in the entire stream.
    std::stringstream ss("\n");
    Buffer b(ss, 128);

    BOOST_REQUIRE_EQUAL(b.getColumn(), 1u);

    while (b.current < b.end || b.load()) {
        ++b.current;
    }
    BOOST_REQUIRE_EQUAL(b.getColumn(), 1u);
}

BOOST_AUTO_TEST_CASE(NewlineCase1)
{
    // Test column position around newlines. This case has a newline as the
    // first char in a new buffer.
    std::stringstream ss("123456\n890\n1");
    Buffer b(ss, 10);
    const int stream_len = 12;
    unsigned int column[stream_len] = {1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 1};

    int ix = 0;
    while (b.current < b.end || b.load()) {
        BOOST_REQUIRE_EQUAL(b.getColumn(), column[ix]);
        ++b.current;
        ++ix;
    }
}

BOOST_AUTO_TEST_SUITE_END()
