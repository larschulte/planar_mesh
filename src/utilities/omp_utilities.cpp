#include "utilities/omp_utilities.hpp"

bool enable_logging = false;

void omp_set_lock_with_log(omp_lock_t &lock, const char *lock_name)
{
    int max_attempts = 10;
    int attempts = 0;

    // First, try using omp_test_lock to see if we can acquire the lock immediately
    while (!omp_test_lock(&lock))
    {
        attempts++;

        // Sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // If max_attempts is exceeded, we can log that we're forcing a wait with omp_set_lock()
        if (attempts >= max_attempts)
        {
            if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " (XXX) " << lock_name << std::endl;
            omp_set_lock(&lock); // Block until the lock is acquired
            return;
        }
    }

    // If we get here, omp_test_lock was successful
    if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " <---> " << lock_name << std::endl;
}

void omp_set_nested_lock_with_log(omp_nest_lock_t &lock, const char *lock_name)
{
    int max_attempts = 10;
    int attempts = 0;

    // First, try using omp_test_nest_lock to see if we can acquire the lock immediately
    while (!omp_test_nest_lock(&lock))
    {
        attempts++;

        // Sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // If max_attempts is exceeded, we can log that we're forcing a wait with omp_set_nest_lock()
        if (attempts >= max_attempts)
        {
            if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " (XXX) " << lock_name << std::endl;
            omp_set_nest_lock(&lock); // Block until the lock is acquired
            return;
        }
    }

    // If we get here, omp_test_nest_lock was successful
    if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " <---> " << lock_name << std::endl;
}

void omp_unset_lock_with_log(omp_lock_t &lock, const char *lock_name)
{
    // Release the lock
    omp_unset_lock(&lock);

    // Log that the lock has been released
    if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << "       " << lock_name << std::endl;
}

void omp_unset_nested_lock_with_log(omp_nest_lock_t &lock, const char *lock_name)
{
    // Release the lock
    omp_unset_nest_lock(&lock);

    // Log that the lock has been released
    if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << "       " << lock_name << std::endl;
}

bool omp_test_lock_with_log(omp_lock_t &lock, const char *lock_name)
{
    // First, try using omp_test_lock to see if we can acquire the lock immediately
    if (!omp_test_lock(&lock))
    {
        if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " test failed '" << lock_name << "'." << std::endl;
        return false;
    }
    else
    {
        if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " test successful '" << lock_name << "'." << std::endl;
        return true;
    }
}

bool omp_test_nested_lock_with_log(omp_nest_lock_t &lock, const char *lock_name)
{
    // First, try using omp_test_nest_lock to see if we can acquire the lock immediately
    if (!omp_test_nest_lock(&lock))
    {
        if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " test failed '" << lock_name << "'." << std::endl;
        return false;
    }
    else
    {
        if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " test successful '" << lock_name << "'." << std::endl;
        return true;
    }
}

// Constructor implementation
OmpLockGuard::OmpLockGuard(omp_lock_t &lock, std::string lock_name)
    : lock_(lock), lock_name_(lock_name)
{
    omp_set_lock_with_log(lock_, lock_name_.c_str());
}

// Destructor implementation
OmpLockGuard::~OmpLockGuard()
{
    omp_unset_lock_with_log(lock_, lock_name_.c_str());
}