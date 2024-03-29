/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// Adapted from Scikit-FMM: https://github.com/scikit-fmm/scikit-fmm

#include "Debug.h"
#include "Heap.h"

/*! \file Heap.cpp
    \brief An implementation of a heap data structure (binary tree).
 */

namespace slsm
{
    Heap::Heap(unsigned int maxLength_, bool isTest_) :
        maxLength(maxLength_),
        isTest(isTest_)
    {
        // Initialise empty heap.
        heapLength = 0;
        listLength = 0;

        // Resize data structures.
        distance.resize(maxLength);
        heap.resize(maxLength);
        address.resize(maxLength);
        backPointer.resize(maxLength);
    }

    unsigned int Heap::push(unsigned int address_, double value)
    {
        // Make sure the heap isn't full.
        errno = 0;
        slsm_check(heapLength < maxLength, "push: Heap is full!");

        // Add entry to the heap.
        heap[heapLength]        = listLength;
        address[listLength]     = address_;
        distance[listLength]    = value;
        backPointer[listLength] = heapLength;

        // Increment heap size.
        heapLength++;

        // Increment list size.
        listLength++;

        // Sift entry down the heap.
        siftDown(0, heapLength-1);

        // Test heap.
        if (isTest) test();

        return listLength - 1;

    error:
        exit(EXIT_FAILURE);
    }

    void Heap::pop(unsigned int& address_, double& value)
    {
        // Make sure the heap isn't empty.
        errno = 0;
        slsm_check(heapLength != 0, "pop: Heap is empty!");

        // Remove entry from heap.
        address_                = address[heap[0]];
        value                   = distance[heap[0]];
        heap[0]                 = heap[heapLength-1];
        backPointer[heap[0]]    = 0;

        // Decrement heap size.
        heapLength--;

        // Sift new heap head up.
        siftUp(0);

        // Test heap.
        if (isTest) test();

        return;

    error:
        exit(EXIT_FAILURE);
    }

    void Heap::siftDown(unsigned int startPos, unsigned int pos)
    {
        // The current item on the heap.
        unsigned int newItem = heap[pos];

        // The index of the parent in the finite element grid.
        unsigned int parent;

        // The position of the parent within the heap.
        unsigned int parentPos;

        while (pos > startPos)
        {
            parentPos = (pos-1) >> 1;
            parent = heap[parentPos];

            if (distance[newItem] < distance[parent])
            {
                // Move item down the heap.
                heap[pos] = parent;
                backPointer[parent] = pos;
                pos = parentPos;

                continue;
            }
            break;
        }

        // Update the heap.
        heap[pos] = newItem;
        backPointer[newItem] = pos;
    }

    void Heap::siftUp(unsigned int pos)
    {
        // Store starting position.
        unsigned int startPos = pos;

        // The current item on the heap.
        unsigned int newItem = heap[pos];

        // Index of right-hand child node.
        unsigned int rightPos;

        // Index of the child node with the minimum distance.
        // Default to the left-hand child.
        unsigned int childPos = 2*pos + 1;

        while (childPos < heapLength)
        {
            rightPos = childPos + 1;

            if (rightPos < heapLength)
            {
                // If right-hand leaf has smaller distance, update the minimum child.
                if (distance[heap[rightPos]] < distance[heap[childPos]])
                    childPos = rightPos;
            }

            // Update the heap.
            heap[pos] = heap[childPos];
            backPointer[heap[childPos]] = pos;
            pos = childPos;

            // Move along the heap to next left-hand child node.
            childPos = 2*pos + 1;
        }

        // Update the heap.
        heap[pos] = newItem;

        // Sift down the heap.
        siftDown(startPos, pos);
    }

    void Heap::set(unsigned int index, double newDistance)
    {
        // Save current distance.
        double oldDistance = distance[index];

        // Update distance.
        distance[index] = newDistance;

        // Current position.
        unsigned int pos = backPointer[index];

        if (newDistance > oldDistance)
        {
            siftUp(pos);
        }

        if (distance[heap[pos]] != newDistance)
        {
            if (isTest) test();
            return;
        }

        siftDown(0, pos);

        if (isTest) test();
    }

    bool Heap::empty() const
    {
        return (heapLength == 0) ? true : false;
    }

    const double& Heap::peek() const
    {
        // Make sure the heap isn't empty.
        errno = 0;
        slsm_check(heapLength != 0, "peek: Heap is empty!");

        return distance[heap[0]];

    error:
        exit(EXIT_FAILURE);
    }

    const unsigned int& Heap::size() const
    {
        return heapLength;
    }

    void Heap::test() const
    {
        // Test heap invariant.
        for (unsigned int i=0;i<heapLength;i++)
        {
            unsigned int children[2];

            children[0] = 2*i + 1;  // Left-hand child.
            children[1] = 2*i + 2;  // Right-hand child.

            for (unsigned int j=0;j<2;j++)
            {
                if (children[j] < (heapLength - 1))
                {
                    double parentDistance = distance[heap[i]];
                    double childDistance = distance[heap[children[j]]];

                    errno = 0;
                    slsm_check(parentDistance <= childDistance, "Heap invariant violation.");
                }
            }
        }

        // Test for backpointer consistency.
        for (unsigned int i=0;i<heapLength;i++)
        {
            slsm_check(backPointer[heap[i]] == i, "Heap backpointer inconsistency.");
        }

        return;

    error:
        exit(EXIT_FAILURE);
    }
}
