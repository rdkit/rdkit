#include "Writer.hpp"

#include <boost/algorithm/string/predicate.hpp>
#ifdef MAEPARSER_HAVE_BOOST_IOSTREAMS
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include "MaeBlock.hpp"

using namespace std;
using boost::algorithm::ends_with;
#ifdef MAEPARSER_HAVE_BOOST_IOSTREAMS
using boost::iostreams::file_sink;
using boost::iostreams::filtering_ostream;
#endif

namespace schrodinger
{
namespace mae
{

Writer::Writer(std::shared_ptr<ostream> stream) : m_out(std::move(stream))
{
    write_opening_block();
}

Writer::Writer(const std::string& fname)
{
    const auto ios_mode = std::ios_base::out | std::ios_base::binary;

    if (ends_with(fname, ".maegz") || ends_with(fname, ".mae.gz")) {
#ifdef MAEPARSER_HAVE_BOOST_IOSTREAMS
        auto* gzip_stream = new filtering_ostream();
        gzip_stream->push(boost::iostreams::gzip_compressor());
        gzip_stream->push(file_sink(fname, ios_mode));
        m_out.reset(static_cast<ostream*>(gzip_stream));
#else
        std::stringstream ss;
        ss << "Unable to open " << fname << " for writing, "
            << "as maeparser was compiled without boost::iostreams support";
        throw std::runtime_error(ss.str());
#endif
    } else {
        auto* file_stream = new ofstream(fname, ios_mode);
        m_out.reset(static_cast<ostream*>(file_stream));
    }

    if (m_out->fail()) {
        std::stringstream ss;
        ss << "Failed to open file \"" << fname << "\" for writing operation.";
        throw std::runtime_error(ss.str());
    }

    write_opening_block();
}

void Writer::write(const std::shared_ptr<Block>& block)
{
    block->write(*m_out);
}

void Writer::write_opening_block()
{
    shared_ptr<Block> b = make_shared<Block>("");
    b->setStringProperty("s_m_m2io_version", "2.0.0");
    write(b);
}

} // namespace mae
} // namespace schrodinger
