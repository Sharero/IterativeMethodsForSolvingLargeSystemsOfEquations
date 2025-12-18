#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
   public:
    using clock = std::chrono::steady_clock;

    Timer() {
        start();
    }

    void start() {
        startTime = clock::now();
        isRunning = true;
    }

    void stop() {
        endTime = clock::now();
        isRunning = false;
    }

    [[nodiscard]]
    double elapsedSeconds() const {
        return std::chrono::duration<double>(
                   (isRunning ? clock::now() : endTime) - startTime)
            .count();
    }

    [[nodiscard]]
    std::int64_t elapsedMilliseconds() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                   (isRunning ? clock::now() : endTime) - startTime)
            .count();
    }

   private:
    clock::time_point startTime{};
    clock::time_point endTime{};

    bool isRunning = false;
};

#endif
