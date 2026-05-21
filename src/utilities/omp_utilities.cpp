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

bool omp_test_nested_lock_with_log(omp_nest_lock_t &lock, const char *lock_name, unsigned int max_attempts)
{
    unsigned int attempts = 0;

    // First, try using omp_test_lock to see if we can acquire the lock immediately
    while (!omp_test_nest_lock(&lock))
    {
        attempts++;

        // Sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::cout << "Thread " << omp_get_thread_num() << " sleep 100 ms for " << lock_name << " with attempts = " << attempts << std::endl;

        // If max_attempts is exceeded, we can log that we're forcing a wait with omp_set_lock()
        if (attempts >= max_attempts)
        {
            if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " (XXX) " << lock_name << std::endl;
            return false;
        }
    }

    // If we get here, omp_test_lock was successful
    if (enable_logging) std::cout << "Thread " << omp_get_thread_num() << " <---> " << lock_name << std::endl;

    return true;
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


namespace custom
{
    custom_lock::custom_lock() : writer_thread_id(), write_lock_count(0) {}

    bool custom_lock::write_locked() const 
    {
        return writer_thread_id != std::thread::id(); 
    }

    bool custom_lock::locked_by_current_thread() const 
    { 
        return writer_thread_id == std::this_thread::get_id(); 
    }

    bool custom_lock::become_write_lock_candidate()
    {
        std::lock_guard<std::mutex> lock_guard(mtx);

        // locked 
        if (write_locked())
        {
            // if locked by current thread, return true
            if (locked_by_current_thread())
            {
                return true;
            }
            // if locked by other thread, return false
            else
            {
                return false;
            }
        }
        else // if not locked
        {
            const bool no_candidate = write_lock_candidate_id == std::thread::id();
            if (no_candidate)
            {
                write_lock_candidate_id = std::this_thread::get_id();
                return true;
            }
            else
            {
                const bool candidate_is_current_thread = write_lock_candidate_id == std::this_thread::get_id();
                if (candidate_is_current_thread)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
    }

    bool custom_lock::set_read_lock()
    {
        std::lock_guard<std::mutex> lock_guard(mtx);

        if (!write_locked())
        {
            // lock shared
            shared_mtx.lock_shared();
            return true;
        }
        else
        {
            if (locked_by_current_thread())
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    void custom_lock::unset_read_lock()
    {    
        // std::lock_guard<std::mutex> lock_guard(mtx);

        if (!write_locked())
        {
            // lock shared
            shared_mtx.unlock_shared();
            return;
        }
        else
        {
            if (locked_by_current_thread())
            {
                return;
            }
            else
            {
                throw std::runtime_error("Thread tried to release a read lock it doesn't own!");
            }
        }
    }

    bool custom_lock::set_write_lock()
    {
        std::lock_guard<std::mutex> lock_guard(mtx);

        if (!write_locked())
        {
            shared_mtx.lock();
            writer_thread_id = std::this_thread::get_id();
            write_lock_count = 1;
            return true;
        }
        else
        {
            if (locked_by_current_thread())
            {
                write_lock_count++;
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    void custom_lock::unset_write_lock()
    {
        std::lock_guard<std::mutex> lock_guard(mtx);

        if (!write_locked())
        {
            throw std::runtime_error("Thread tried to release a write lock that is not locked");
        }
        else
        {
            if (locked_by_current_thread())
            {
                write_lock_count--;
                if (write_lock_count == 0)
                {
                    writer_thread_id = std::thread::id();
                    write_lock_candidate_id = std::thread::id();
                    shared_mtx.unlock();
                }
            }
            else
            {
                std::cout << "Thread " << std::this_thread::get_id() << " tried to release a write lock it doesn't own!" << std::endl;
                std::cout << "Thread " << writer_thread_id << " is the owner of the write lock." << std::endl;
                std::cout << "Thread " << write_lock_candidate_id << " is the write lock candidate." << std::endl;
                throw std::runtime_error("Thread tried to release a write lock it doesn't own!");
            }
        }
    }
} // namespace custom