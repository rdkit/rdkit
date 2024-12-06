// TEMPORARY PLACE HOLDER
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

    // Returns the current Unicode code point without advancing the iterator
    uint32_t GetCharacter() const {
        if (AtEnd())
            throw std::out_of_range("Iterator has reached the end of the string.");
        return currentCodePoint;
    }

    // Checks if the iterator has reached the end of the string
    bool AtEnd() const {
        return pos >= data.size();
    }

    // Post-increment operator: advances to the next character
    UTF8Iterator operator++(int) {
        UTF8Iterator temp = *this; // Copy current state
        if (!AtEnd())
            currentCodePoint = DecodeNext();
        return temp; // Return old state
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
