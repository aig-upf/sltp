/*
  Copyright (C) 2011 by the PSVN Research Group, University of Alberta
*/

#ifndef _PRIORITY_QUEUE_HPP
#define _PRIORITY_QUEUE_HPP

#include <queue>
#include <cstdio>
#include <cassert>
#include <cstdlib>

template <typename T>
class PriorityQueue
{
public:
    PriorityQueue();
    int Add(int f, int g, const T& data);
    const T& Top() const;
    int CurrentPriority() const;
    bool Empty() const;
    void Pop();
    void Clear();

    void Modify(int f, int g, int index, const T& data);

private:
    /* Min priority queue on f. 
       If same f: prioritize on max g. */
    class Bucket
    {
    public:
        void Pop();
        bool Empty() const;
        int Add(int g, const T& data);
        const T& Top() const;
        void Clear();
        void Modify(int g, int index, const T& data);
    private:
        std::vector<std::vector<T> > m_stacks;
    };
    int m_max;
    int m_at;
    std::vector<Bucket> m_buckets;
};

template<typename T>
const T& PriorityQueue<T>::Bucket::Top() const
{
    assert(!m_stacks.empty());
    return m_stacks.back().back();
}

template<typename T>
void PriorityQueue<T>::Bucket::Pop()
{
    m_stacks.back().pop_back();
    while (!m_stacks.empty() && m_stacks.back().empty())
        m_stacks.pop_back();
}

template<typename T>
bool PriorityQueue<T>::Bucket::Empty() const
{
    return m_stacks.empty();
}

template<typename T>
int PriorityQueue<T>::Bucket::Add(int g, const T& data)
{
    if (m_stacks.size() <= (size_t)g)
        m_stacks.resize(g + 1);
    m_stacks[g].push_back(data);
    return m_stacks[g].size() - 1;
}

template<typename T>
void PriorityQueue<T>::Bucket::Modify(int g, int index, const T& data)
{
    // if (g >= (int)m_stacks.size() || index >= m_stacks[g].size()) {
    //     printf("%d %d | %d %d\n", g, index, m_stacks.size(), m_stacks[g].size());
    // }
    assert(g < (int)m_stacks.size());
    assert(index < (int)m_stacks[g].size());
    m_stacks[g][index] = data;
}

template<typename T>
void PriorityQueue<T>::Bucket::Clear()
{
    return m_stacks.clear();
}

template<typename T>
PriorityQueue<T>::PriorityQueue()
    : m_max(128),
      m_at(m_max),
      m_buckets(m_max)
{
}

template<typename T>
void PriorityQueue<T>::Modify(int f, int g, int index, const T& data)
{
    assert(f < (int)m_buckets.size());
    m_buckets[f].Modify(g, index, data);
}

template<typename T>
int PriorityQueue<T>::Add(int f, int g, const T& data)
{
    assert(f >= 0);
    if (f >= m_max) {
        // double size until we can accomodate f
        int new_max = m_max;
        while (f >= new_max)
            new_max *= 2;
        if (m_at == m_max)
            m_at = new_max;
        m_max = new_max;
        m_buckets.resize(new_max);
    }
    if (f < m_at)
        m_at = f;
    return m_buckets[f].Add(g, data);
}

template<typename T>
const T& PriorityQueue<T>::Top() const
{
    return m_buckets[m_at].Top();
}

template<typename T>
int PriorityQueue<T>::CurrentPriority() const
{
    return m_at;
}

template<typename T>
bool PriorityQueue<T>::Empty() const
{
    return m_at == m_max;
}

template<typename T>
void PriorityQueue<T>::Pop()
{
    m_buckets[m_at].Pop();
    while (m_at < m_max && m_buckets[m_at].Empty())
        ++m_at;
}

template<typename T>
void PriorityQueue<T>::Clear()
{
    for (int i = 0; i < m_max; ++i)
        m_buckets[i].Clear();
    m_at = m_max;
}

#endif // _PRIORITY_QUEUE_HPP
