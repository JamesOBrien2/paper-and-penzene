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
// Modified for Penzene: the original returned the whole string's length from
// GetUTF8Length() and the index after the character from GetByteIndex(), so
// MakeStringSafe copied a growing suffix per character (garbled CDX/CDXML text).
// This version tracks the current character's start and byte length, stays within
// the string, and reads a byte that doesn't start valid UTF-8 as a one-byte character
// (what MakeStringSafe expects, to map legacy single-byte symbols).
#ifndef UTF8_ITERATOR_H
#define UTF8_ITERATOR_H

#include <cstdint>
#include <stdexcept>
#include <string>

class UTF8Iterator {
public:
    explicit UTF8Iterator(const std::string& str) : data(str) { Decode(); }

    size_t GetUTF8Length() const { return length; }   // bytes in the current character
    size_t GetByteIndex() const { return start; }     // where it starts
    bool AtEnd() const { return start >= data.size(); }

    // The current Unicode code point
    uint32_t GetCharacter() const {
        if (AtEnd()) throw std::out_of_range("Iterator has reached the end of the string.");
        return codePoint;
    }

    UTF8Iterator& operator++() {
        if (!AtEnd()) start += length, Decode();
        return *this;
    }
    UTF8Iterator operator++(int) {
        UTF8Iterator old = *this;
        ++*this;
        return old;
    }
    void Next() { ++*this; }

private:
    const std::string& data;
    size_t start = 0, length = 0;
    uint32_t codePoint = 0;

    void Decode() {
        if (AtEnd()) {
            length = 0, codePoint = 0;
            return;
        }
        const unsigned char lead = data[start];
        const size_t need = lead < 0x80 ? 1 : (lead & 0xE0) == 0xC0 ? 2 : (lead & 0xF0) == 0xE0 ? 3 : (lead & 0xF8) == 0xF0 ? 4 : 0;
        bool valid = need > 0 && start + need <= data.size();
        for (size_t k = 1; valid && k < need; ++k) valid = (static_cast<unsigned char>(data[start + k]) & 0xC0) == 0x80;
        if (!valid) {  // not UTF-8 here: one byte, its own value
            length = 1, codePoint = lead;
            return;
        }
        length = need;
        codePoint = need == 1 ? lead : need == 2 ? lead & 0x1F : need == 3 ? lead & 0x0F : lead & 0x07;
        for (size_t k = 1; k < need; ++k) codePoint = (codePoint << 6) | (static_cast<unsigned char>(data[start + k]) & 0x3F);
    }
};
#endif
