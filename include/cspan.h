#ifndef CSPAN_H
#define CSPAN_H

#include <vector>

#include "../include/span.h"

class CSpan {
   public:
    const std::vector<double>& vec;
    const size_t start;
    const size_t size;

    CSpan(const std::vector<double>& vec, size_t start, size_t size)
        : vec{vec}, start{start}, size{size} {}

    inline CSpan(const Span& span)
        : vec{span.vec}, start{span.start}, size{span.size} {}

    inline auto begin() const -> const double* {
        return vec.data() + start;
    }
};

#endif