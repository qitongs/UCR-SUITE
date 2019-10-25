//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_CIRCULAR_ARRAY_H
#define UCR_SUITE_CIRCULAR_ARRAY_H

/// TODO depreciate this file by boost's circular buffer

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct CircularArray
{
    int * dq;
    int size, capacity;
    int f, r;
};

/// Initial the queue at the beginning step of envelop calculation
void init(CircularArray * d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int) *d->capacity);
    d->f = 0;
    d->r = d->capacity-1;
}

/// Destroy the queue
void destroy(CircularArray * d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct CircularArray * d, int v)
{
    d->dq[d->r] = v;
    d->r --;
    if (d->r < 0)
        d->r = d->capacity - 1;
    d->size ++;
}

/// Delete the current (front) element from queue
void pop_front(struct CircularArray * d)
{
    d->f --;
    if (d->f < 0)
        d->f = d->capacity - 1;
    d->size --;
}

/// Delete the last element from queue
void pop_back(struct CircularArray * d)
{
    d->r = (d->r + 1) % d->capacity;
    d->size --;
}

/// Get the value at the current position of the circular queue
int front(struct CircularArray * d)
{
    int aux = d->f - 1;
    if (aux < 0)
        aux = d->capacity - 1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct circular_array *d)
int back(struct CircularArray * d)
{
    int aux = (d->r + 1) % d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct CircularArray * d)
{
    return d->size == 0;
}

#endif //UCR_SUITE_CIRCULAR_ARRAY_H
