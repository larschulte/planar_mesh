#pragma once

#include <omp.h>
#include <iostream>
#include <thread>  // for std::this_thread::sleep_for
#include <chrono>  // for std::chrono::milliseconds

void omp_set_lock_with_log(omp_lock_t &lock, const char* lock_name);
void omp_set_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);

void omp_unset_lock_with_log(omp_lock_t &lock, const char* lock_name);
void omp_unset_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);

bool omp_test_lock_with_log(omp_lock_t &lock, const char* lock_name);
bool omp_test_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name);