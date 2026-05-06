#pragma once

#include <string>
#include <sstream>

#include <boost/test/unit_test.hpp>

template <class T> class CheckExceptionMsg
{
  public:
    CheckExceptionMsg(std::string msg)
        : m_expected(std::move(msg)), m_match(false)
    {
    }

    ~CheckExceptionMsg()
    {
        // This message should be printer after BOOST_REQUIRE_EXCEPTION's
        std::stringstream ss;
        ss << '\"' << m_expected << "\" not found in \"" << m_message << '\"';

        BOOST_CHECK_MESSAGE(m_match, ss.str());
    }

    bool operator()(const T& exc)
    {
        m_message = exc.what();
        m_match = m_message.find(m_expected) != std::string::npos;

        return m_match;
    }

  private:
    const std::string m_expected;
    std::string m_message;
    bool m_match;
};