#pragma once

#include <stack>
#include <queue>
#include <memory>

// Common interface
template <typename T>
class DataStructure
{
public:
    virtual void push(const T& value) = 0;
    virtual void pop() = 0;
    virtual T get() = 0;
    virtual bool empty() const = 0;
    virtual unsigned int size() const = 0;
    virtual ~DataStructure() = default;
};

// Concrete class for Queue
template <typename T>
class QueueStructure : public DataStructure<T>
{
private:
    std::queue<T> queue;

public:
    void push(const T& value) override;
    void pop() override;
    T get() override;
    bool empty() const override;
    unsigned int size() const override;
};

// Concrete class for Stack
template <typename T>
class StackStructure : public DataStructure<T>
{
private:
    std::stack<T> stack;

public:
    void push(const T& value) override;
    void pop() override;
    T get() override;
    bool empty() const override;
    unsigned int size() const override;
};

// Wrapper class that decides to use stack or queue
template <typename T>
class queue_or_stack
{
private:
    std::unique_ptr<DataStructure<T>> structure;

public:
    // Constructor decides whether to use queue or stack
    queue_or_stack() : structure(nullptr) {};
    queue_or_stack(bool useQueue);

    void push(const T& value);
    void pop();
    T get();
    bool empty() const;
    unsigned int size() const;
};

#include "utilities/queue_or_stack.tpp"