#ifndef SPAN_H
#define SPAN_H

#include <vector>

class Span {
   public:
    std::vector<double>& vec;
    const size_t start;
    const size_t size;

    Span(std::vector<double>& vec, size_t start, size_t size)
        : vec{vec}, start{start}, size{size} {}

    inline auto begin() const -> double* {
        return vec.data() + start;
    }
};

#endif