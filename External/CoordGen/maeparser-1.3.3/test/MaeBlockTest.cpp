#include <iostream>
#include <memory>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include "MaeBlock.hpp"
#include "MaeConstants.hpp"

using namespace schrodinger;

BOOST_AUTO_TEST_SUITE(MaeBlockSuite)

BOOST_AUTO_TEST_CASE(maeBlock)
{
    double tolerance = std::numeric_limits<double>::epsilon();

    mae::Block b("dummy");
    b.setRealProperty("a", 1.0);
    BOOST_REQUIRE(b.hasRealProperty("a"));
    BOOST_REQUIRE(!b.hasRealProperty("b"));
    BOOST_REQUIRE_CLOSE(b.getRealProperty("a"), 1.0, tolerance);
    BOOST_REQUIRE_THROW(b.getRealProperty("b"), std::out_of_range);

    b.setIntProperty("a", 3);
    BOOST_REQUIRE(b.hasIntProperty("a"));
    BOOST_REQUIRE(!b.hasIntProperty("b"));
    BOOST_REQUIRE_EQUAL(b.getIntProperty("a"), 3);
    BOOST_REQUIRE_THROW(b.getIntProperty("b"), std::out_of_range);

    b.setBoolProperty("a", true);
    BOOST_REQUIRE(b.hasBoolProperty("a"));
    BOOST_REQUIRE(!b.hasBoolProperty("b"));
    BOOST_REQUIRE(b.getBoolProperty("a"));
    BOOST_REQUIRE_THROW(b.getBoolProperty("b"), std::out_of_range);

    const std::vector<std::string> strings = {"Regular", "Spaced String"};
    for (const auto& value : strings) {
        b.setStringProperty("a", value);
        BOOST_REQUIRE(b.hasStringProperty("a"));
        BOOST_REQUIRE(!b.hasStringProperty("b"));
        BOOST_REQUIRE_EQUAL(b.getStringProperty("a"), value);
        BOOST_REQUIRE_THROW(b.getStringProperty("b"), std::out_of_range);
    }

    // All the properties we set have the same name, but different types
    const std::string property_name("a");

    auto bool_props = b.getProperties<mae::BoolProperty>();
    BOOST_CHECK_EQUAL(bool_props.size(), 1u);
    BOOST_CHECK_EQUAL(bool_props.begin()->first, property_name);

    auto int_props = b.getProperties<int>();
    BOOST_CHECK_EQUAL(int_props.size(), 1u);
    BOOST_CHECK_EQUAL(int_props.begin()->first, property_name);

    auto double_props = b.getProperties<double>();
    BOOST_CHECK_EQUAL(double_props.size(), 1u);
    BOOST_CHECK_EQUAL(double_props.begin()->first, property_name);

    auto string_props = b.getProperties<std::string>();
    BOOST_CHECK_EQUAL(string_props.size(), 1u);
    BOOST_CHECK_EQUAL(string_props.begin()->first, property_name);
}

BOOST_AUTO_TEST_CASE(maeIndexedRealProperty)
{
    double tolerance = std::numeric_limits<double>::epsilon();
    {
        auto dv = std::make_shared<std::vector<double>>();
        dv->push_back(1.0);
        dv->push_back(2.0);
        dv->push_back(3.0);
        mae::IndexedRealProperty irp(*dv);
        BOOST_REQUIRE_CLOSE(irp[0], 1.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp[1], 2.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp[2], 3.0, tolerance);
    }
    {
        std::vector<double> dv;
        dv.push_back(1.0);
        dv.push_back(2.0);
        dv.push_back(3.0);
        mae::IndexedRealProperty irp(dv);
        BOOST_REQUIRE_CLOSE(irp[0], 1.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp[1], 2.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp[2], 3.0, tolerance);
    }
    {
        auto dv = std::make_shared<std::vector<double>>();
        boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
        bs->set(1);

        dv->push_back(1.0);
        dv->push_back(0.0);
        dv->push_back(3.0);
        mae::IndexedRealProperty irp(*dv, bs);
        BOOST_REQUIRE(irp.isDefined(0));
        BOOST_REQUIRE_CLOSE(irp[0], 1.0, tolerance);
        BOOST_REQUIRE(!irp.isDefined(1));
        BOOST_REQUIRE_THROW(irp[1], std::runtime_error);
        BOOST_REQUIRE(irp.isDefined(2));
        BOOST_REQUIRE_CLOSE(irp[2], 3.0, tolerance);
    }
    {
        std::vector<double> dv;
        boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
        bs->set(1);

        dv.push_back(1.0);
        dv.push_back(0.0);
        dv.push_back(3.0);
        mae::IndexedRealProperty irp(dv, bs);
        BOOST_REQUIRE(irp.isDefined(0));
        BOOST_REQUIRE_CLOSE(irp[0], 1.0, tolerance);
        BOOST_REQUIRE(!irp.isDefined(1));
        BOOST_REQUIRE_THROW(irp[1], std::runtime_error);
        BOOST_REQUIRE(irp.isDefined(2));
        BOOST_REQUIRE_CLOSE(irp[2], 3.0, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(maeIndexedBlock)
{
    using namespace mae;
    double tolerance = std::numeric_limits<double>::epsilon();
    {
        std::vector<double> dv;
        boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
        bs->set(1);

        dv.push_back(1.0);
        dv.push_back(0.0);
        dv.push_back(3.0);
        IndexedBlock ib("m_atom");
        BOOST_REQUIRE(!ib.hasRealProperty("r_m_float"));
        auto irps = std::make_shared<IndexedRealProperty>(dv, bs);

        ib.setRealProperty("r_m_float", irps);
        BOOST_REQUIRE(ib.hasRealProperty("r_m_float"));

        auto irpg = ib.getRealProperty("r_m_float");
        IndexedRealProperty& irp = *irpg;
        BOOST_REQUIRE(irp.isDefined(0));
        BOOST_REQUIRE_CLOSE(irp[0], 1.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp.at(0, 999.0), 1.0, tolerance);
        BOOST_REQUIRE(!irp.isDefined(1));
        BOOST_REQUIRE_THROW(irp[1], std::runtime_error);
        BOOST_REQUIRE_CLOSE(irp.at(1, 999.0), 999.0, tolerance);
        BOOST_REQUIRE(irp.isDefined(2));
        BOOST_REQUIRE_CLOSE(irp[2], 3.0, tolerance);
        BOOST_REQUIRE_CLOSE(irp.at(2, 999.0), 3.0, tolerance);

        auto no_prop = ib.getRealProperty("r_m_nonexistant");
        BOOST_REQUIRE(no_prop.get() == nullptr);
    }
}

BOOST_AUTO_TEST_CASE(maeIndexedBlockBool)
{
    using namespace mae;
    {
        std::vector<BoolProperty> dv;
        boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
        bs->set(1);

        dv.push_back(true);
        dv.push_back(false);
        dv.push_back(true);
        IndexedBlock ib("m_atom");
        BOOST_REQUIRE(!ib.hasBoolProperty("b_m_bool"));
        auto ibps = std::make_shared<IndexedBoolProperty>(dv, bs);

        ib.setBoolProperty("b_m_bool", ibps);
        BOOST_REQUIRE(ib.hasBoolProperty("b_m_bool"));

        auto ibpg = ib.getBoolProperty("b_m_bool");
        IndexedBoolProperty& ibp = *ibpg;
        BOOST_REQUIRE(ibp.isDefined(0));
        BOOST_REQUIRE_EQUAL(ibp[0], static_cast<BoolProperty>(true));
        BOOST_REQUIRE(!ibp.isDefined(1));
        BOOST_REQUIRE_THROW(ibp[1], std::runtime_error);
        BOOST_REQUIRE(ibp.isDefined(2));
        BOOST_REQUIRE_EQUAL(ibp[2], static_cast<BoolProperty>(true));
    }
}

BOOST_AUTO_TEST_CASE(maeIndexedBlockString)
{
    using namespace mae;
    {
        std::vector<std::string> dv;
        boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
        bs->set(1);

        dv.push_back("Hi with space");
        dv.push_back("ignore me");
        dv.push_back("Bye");
        IndexedBlock ib("m_atom");
        BOOST_REQUIRE(!ib.hasStringProperty("s_m_string"));
        auto isps = std::make_shared<IndexedStringProperty>(dv, bs);

        ib.setStringProperty("s_m_string", isps);
        BOOST_REQUIRE(ib.hasStringProperty("s_m_string"));

        auto ispg = ib.getStringProperty("s_m_string");
        IndexedStringProperty& isp = *ispg;
        BOOST_REQUIRE(isp.isDefined(0));
        BOOST_REQUIRE_EQUAL(isp[0], "Hi with space");
        BOOST_REQUIRE(!isp.isDefined(1));
        BOOST_REQUIRE_THROW(isp[1], std::runtime_error);
        BOOST_REQUIRE(isp.isDefined(2));
        BOOST_REQUIRE_EQUAL(isp[2], "Bye");
    }
}

std::shared_ptr<mae::IndexedBlock> getExampleIndexedBlock()
{
    using namespace mae;
    auto ib = std::make_shared<IndexedBlock>("m_atom");

    // Set up a bool property
    std::vector<BoolProperty> dv = {true, false, true};
    boost::dynamic_bitset<>* bs = new boost::dynamic_bitset<>(3);
    bs->set(1);

    auto ibps = std::make_shared<IndexedBoolProperty>(dv, bs);
    ib->setBoolProperty("b_m_bool", ibps);

    // Set up a real property
    std::vector<double> rv = {0.1, 42};
    boost::dynamic_bitset<>* rbs = new boost::dynamic_bitset<>(3);
    rbs->set(2);

    auto irps = std::make_shared<IndexedRealProperty>(rv, rbs);
    ib->setRealProperty("r_m_reals", irps);

    return ib;
}

BOOST_AUTO_TEST_CASE(toStringProperties)
{
    const std::string rval =
        R"(dummy {
  b_m_bool
  r_m_real
  i_m_int
  s_m_string
  :::
  0
  1.000000
  42
  "mae\"parser"
  m_atom[3] {
    # First column is Index #
    b_m_bool
    r_m_reals
    :::
    1 1 0.100000
    2 <> 42.000000
    3 1 <>
    :::
  }
}

)";
    mae::Block b("dummy");
    b.setRealProperty("r_m_real", 1.0);
    b.setBoolProperty("b_m_bool", false);
    b.setIntProperty("i_m_int", 42);
    b.setStringProperty("s_m_string", "mae\"parser");

    auto ib = getExampleIndexedBlock();

    auto ibm = std::make_shared<mae::IndexedBlockMap>();
    ibm->addIndexedBlock(ib->getName(), ib);
    b.setIndexedBlockMap(ibm);

    BOOST_REQUIRE_EQUAL(b.toString(), rval);

    // All the properties we set have the same name, but different types

    auto bool_props = b.getProperties<mae::BoolProperty>();
    BOOST_CHECK_EQUAL(bool_props.size(), 1u);
    BOOST_CHECK_EQUAL(bool_props.begin()->first, std::string("b_m_bool"));

    auto int_props = b.getProperties<int>();
    BOOST_CHECK_EQUAL(int_props.size(), 1u);
    BOOST_CHECK_EQUAL(int_props.begin()->first, std::string("i_m_int"));

    auto double_props = b.getProperties<double>();
    BOOST_CHECK_EQUAL(double_props.size(), 1u);
    BOOST_CHECK_EQUAL(double_props.begin()->first, std::string("r_m_real"));

    auto string_props = b.getProperties<std::string>();
    BOOST_CHECK_EQUAL(string_props.size(), 1u);
    BOOST_CHECK_EQUAL(string_props.begin()->first, std::string("s_m_string"));
}

BOOST_AUTO_TEST_CASE(toStringIndexedProperties)
{
    using namespace mae;
    const std::string rval =
        R"(m_atom[3] {
  # First column is Index #
  b_m_bool
  r_m_reals
  :::
  1 1 0.100000
  2 <> 42.000000
  3 1 <>
  :::
}
)";

    auto ib = getExampleIndexedBlock();

    BOOST_REQUIRE_EQUAL(ib->toString(), rval);
}

BOOST_AUTO_TEST_SUITE_END()
