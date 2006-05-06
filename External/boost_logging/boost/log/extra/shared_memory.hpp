// shared_memory.hpp

// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 John Torjo (john@torjo.com)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef JT_BOOST_LOG_shared_memory_HPP
#define JT_BOOST_LOG_shared_memory_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif


#include <boost/log/log_impl.hpp>
#include <boost/shmem/shmem_named_shared_object.hpp>
#include <boost/log/detail/ts_none.hpp>
#include <boost/shared_ptr.hpp>


namespace boost { namespace logging {

template<class mutex_type = logging_types::mutex, class lock_type = logging_types::lock> struct basic_shared_memory_appender {
    typedef logging_types::char_t char_t;

    enum { just_in_case = 8192 };
    basic_shared_memory_appender(const std::string & name, std::size_t mem_size = 2 * 1024 * 1024 * sizeof(char_t) ) : m_info(new info) {
        m_info->mem_size = mem_size;
        
        // this segment might have been previously created...
        m_info->segment.open_or_create(name.c_str(), m_info->mem_size + just_in_case);

        // the string
        { typedef std::pair<char_t*, std::size_t> pair;
        pair res = m_info->segment.find<char_t>("shared_log_object");
        if ( !res.first)
            // we're creating it right now
            m_info->segment.construct<char_t>("shared_log_object")[mem_size](0);

        res = m_info->segment.find<char_t>("shared_log_object");
        assert( res.first); // should be created by now
        m_info->memory = res.first;
        }

        // the occupied size
        { typedef std::pair<long*, std::size_t> pair;
        pair res = m_info->segment.find<long>("shared_occupied_size");
        if ( !res.first) 
            // we're creating it right now
            m_info->segment.construct<long>("shared_occupied_size")[1](0);

        res = m_info->segment.find<long>("shared_occupied_size");
        assert( res.first); // should be created by now
        m_info->occupied_size = res.first;
        }
    }
    
    void operator () (const logging_types::string&, const logging_types::string& msg) {
        lock_type lk(m_info->mutex);

        bool can_fit = *m_info->occupied_size + msg.size() < m_info->mem_size;
        if ( can_fit) {
            std::copy(msg.begin(), msg.end(), m_info->memory + *m_info->occupied_size);
            *m_info->occupied_size += (long)msg.size();
        }
        else {
            // exceeds bounds
            if ( msg.size() < m_info->mem_size) {
                // move what was previously written, to the left, to make room
                std::size_t keep = m_info->mem_size / 2;
                if ( keep + msg.size() > m_info->mem_size) keep = m_info->mem_size - msg.size();
                std::copy_backward( 
                    m_info->memory + *m_info->occupied_size - keep,
                    m_info->memory + *m_info->occupied_size,
                    m_info->memory + keep);
                std::copy( msg.begin(), msg.end(), m_info->memory + keep);
                *m_info->occupied_size = (long)(keep + msg.size());
            }
            else {
                // message too big
                std::copy(msg.begin(), msg.begin() + m_info->mem_size, m_info->memory);
                *m_info->occupied_size = (long)m_info->mem_size;
            }
        }
    }

private:
    struct info {
        info() : occupied_size(0), memory(0), mem_size(0) {
            // note: we don't want to destroy this segment, since we want it to outlive our
            // application. In the case of a problem, we need to have another program that will
            // take a look at what we just logged.
        }

        mutex_type mutex;
        // how much of the memory is occupied?
        long * occupied_size;
        // the shared memory
        char_t * memory;
        boost::shmem::named_shared_object segment;
        std::size_t mem_size;
    
    };
    boost::shared_ptr<info> m_info;

};


typedef basic_shared_memory_appender<threading::no_mutex, threading::no_lock> shared_memory_appender_single_thread;
typedef basic_shared_memory_appender<> shared_memory_appender;

}}


#endif

