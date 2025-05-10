// Added by Glydade: Reimplementation of CDMap

// BSD 3-Clause License
// 
// Copyright (c) 2025, Glysade Inc
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef UTF8_ITERATOR_H
#define UTF8_ITERATOR_H

#include <string>
#include <stdexcept>
#include <cstdint>
#include <iterator>
#include <iostream>
#include <locale>
#include <codecvt>


class UTF8Iterator {
public:
    UTF8Iterator(const std::string& str)
      : data(str), pos(0), currentCodePoint(DecodeNext()) {}

    size_t      GetUTF8Length() const {
      UTF8Iterator iter(data);
      int count = 0;
      for(; !iter.AtEnd(); iter++, count++) {}
      return count;
    }

    size_t GetByteIndex() const { return pos; }
  
    // Returns the current Unicode code point without advancing the iterator
    uint32_t GetCharacter() const {
        if (AtEnd())
            throw std::out_of_range("Iterator has reached the end of the string.");
        return currentCodePoint;
    }

    // Checks if the iterator has reached the end of the string
    bool AtEnd() const {
        return pos > data.size();
    }

    // Post-increment operator: advances to the next character
    UTF8Iterator operator++(int) {
        UTF8Iterator temp = *this; // Copy current state
        if (!AtEnd())
            currentCodePoint = DecodeNext();
        return temp; // Return old state
    }
  
    UTF8Iterator& operator++() {
        if (!AtEnd())
            currentCodePoint = DecodeNext();
        return *this;
    }
    

private:
    const std::string& data;
    size_t pos;
    uint32_t currentCodePoint;

    // Decodes the next UTF-8 sequence and advances the position
    uint32_t DecodeNext() {
        if (AtEnd())
            return 0;

        uint32_t codePoint = 0;
        unsigned char byte1 = data[pos];

        if (byte1 < 0x80) { 
            codePoint = byte1;
            ++pos;
        } else if ((byte1 & 0xE0) == 0xC0) { 
            codePoint = (byte1 & 0x1F) << 6;
            codePoint |= (data[pos + 1] & 0x3F);
            pos += 2;
        } else if ((byte1 & 0xF0) == 0xE0) { 
            codePoint = (byte1 & 0x0F) << 12;
            codePoint |= (data[pos + 1] & 0x3F) << 6;
            codePoint |= (data[pos + 2] & 0x3F);
            pos += 3;
        } else if ((byte1 & 0xF8) == 0xF0) { 
            codePoint = (byte1 & 0x07) << 18;
            codePoint |= (data[pos + 1] & 0x3F) << 12;
            codePoint |= (data[pos + 2] & 0x3F) << 6;
            codePoint |= (data[pos + 3] & 0x3F);
            pos += 4;
        } else {
            throw std::runtime_error("Invalid UTF-8 sequence.");
        }

        return codePoint;
    }
};
#endif
