#include "utilities/queue_or_stack.hpp"

//
// QueueStructure
// 

template <typename T>
void QueueStructure<T>::push(T value)
{
    queue.push(value);
}

template <typename T>
void QueueStructure<T>::pop()
{
    if (!queue.empty())
    {
        queue.pop();
    }
}

template <typename T>
T QueueStructure<T>::get()
{
    return queue.front();
}

template <typename T>
bool QueueStructure<T>::empty() const
{
    return queue.empty();
}

template <typename T>
unsigned int QueueStructure<T>::size() const
{
    return queue.size();
}

//
// StackStructure
//

template <typename T>
void StackStructure<T>::push(T value)
{
    stack.push(value);
}

template <typename T>
void StackStructure<T>::pop()
{
    if (!stack.empty())
    {
        stack.pop();
    }
}

template <typename T>
T StackStructure<T>::get()
{
    return stack.top();
}

template <typename T>
bool StackStructure<T>::empty() const
{
    return stack.empty();
}

template <typename T>
unsigned int StackStructure<T>::size() const
{
    return stack.size();
}

// DynamicStructure method definitions
template <typename T>
queue_or_stack<T>::queue_or_stack(bool useQueue)
{
    if (useQueue)
    {
        structure = std::make_unique<QueueStructure<T>>();
    }
    else
    {
        structure = std::make_unique<StackStructure<T>>();
    }
}

template <typename T>
void queue_or_stack<T>::push(T value)
{
    structure->push(value);
}

template <typename T>
void queue_or_stack<T>::pop()
{
    structure->pop();
}

template <typename T>
T queue_or_stack<T>::get()
{
    return structure->get();
}

template <typename T>
bool queue_or_stack<T>::empty() const
{
    return structure->empty();
}

template <typename T>
unsigned int queue_or_stack<T>::size() const
{
    return structure->size();
}