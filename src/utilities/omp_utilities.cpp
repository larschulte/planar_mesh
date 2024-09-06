#include "utilities/omp_utilities.hpp"

void omp_set_lock_with_log(omp_lock_t &lock, const char* lock_name) {
    int max_attempts = 10;
    int attempts = 0;
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // First, try using omp_test_lock to see if we can acquire the lock immediately
    while (!omp_test_lock(&lock)) {
        attempts++;
        
        // Sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // If max_attempts is exceeded, we can log that we're forcing a wait with omp_set_lock()
        if (attempts >= max_attempts) {
            std::cout << "Thread " << thread_id << " (XXX) " << lock_name << std::endl;
            omp_set_lock(&lock);  // Block until the lock is acquired
            return;
        }
    }

    // If we get here, omp_test_lock was successful
    std::cout << "Thread " << thread_id << " <---> " << lock_name << std::endl;
}

void omp_set_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name) {
    int max_attempts = 10;
    int attempts = 0;
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // First, try using omp_test_nest_lock to see if we can acquire the lock immediately
    while (!omp_test_nest_lock(&lock)) {
        attempts++;
        
        // Sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // If max_attempts is exceeded, we can log that we're forcing a wait with omp_set_nest_lock()
        if (attempts >= max_attempts) {
            std::cout << "Thread " << thread_id << " (XXX) " << lock_name << std::endl;
            omp_set_nest_lock(&lock);  // Block until the lock is acquired
            return;
        }
    }

    // If we get here, omp_test_nest_lock was successful
    std::cout << "Thread " << thread_id << " <---> " << lock_name << std::endl;
}

void omp_unset_lock_with_log(omp_lock_t &lock, const char* lock_name) {
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // Release the lock
    omp_unset_lock(&lock);

    // Log that the lock has been released
    std::cout << "Thread " << thread_id << "       " << lock_name << std::endl;
}

void omp_unset_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name) {
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // Release the lock
    omp_unset_nest_lock(&lock);

    // Log that the lock has been released
    std::cout << "Thread " << thread_id << "       " << lock_name << std::endl;
}

bool omp_test_lock_with_log(omp_lock_t &lock, const char* lock_name) {
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // First, try using omp_test_lock to see if we can acquire the lock immediately
    if (!omp_test_lock(&lock)) 
    {
        std::cout << "Thread " << thread_id << " test failed '" << lock_name << "'." << std::endl;
        return false;
    } 
    else 
    {
        std::cout << "Thread " << thread_id << " test successful '" << lock_name << "'." << std::endl;
        return true;
    }
}

bool omp_test_nested_lock_with_log(omp_nest_lock_t &lock, const char* lock_name) {
    int thread_id = omp_get_thread_num();  // Get the thread number inside the function

    // First, try using omp_test_nest_lock to see if we can acquire the lock immediately
    if (!omp_test_nest_lock(&lock)) 
    {
        std::cout << "Thread " << thread_id << " test failed '" << lock_name << "'." << std::endl;
        return false;
    } 
    else 
    {
        std::cout << "Thread " << thread_id << " test successful '" << lock_name << "'." << std::endl;
        return true;
    }
}