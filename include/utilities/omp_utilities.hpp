#pragma once

#include <omp.h>
#include <iostream>
#include <thread>  // for std::this_thread::sleep_for
#include <chrono>  // for std::chrono::milliseconds

// Global flag to control logging at runtime
extern bool enable_logging;

void omp_set_lock_with_log(omp_lock_t &lock, const char* lock_name);
void omp_set_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);

void omp_unset_lock_with_log(omp_lock_t &lock, const char* lock_name);
void omp_unset_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);

bool omp_test_lock_with_log(omp_lock_t &lock, const char* lock_name);
bool omp_test_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);

class OmpLockGuard
{
public:
    // Constructor to initialize and lock the omp_lock_t
    OmpLockGuard(omp_lock_t &lock, std::string lock_name);

    // Destructor to automatically unlock the omp_lock_t when out of scope
    ~OmpLockGuard();

    // Delete copy constructor and assignment operator to prevent copying
    OmpLockGuard(const OmpLockGuard &) = delete;
    OmpLockGuard &operator=(const OmpLockGuard &) = delete;

private:
    omp_lock_t &lock_;      // Reference to the OpenMP lock
    std::string lock_name_; // Name of the lock (for logging)
};