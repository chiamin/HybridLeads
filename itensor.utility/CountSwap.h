#ifndef __COUNT_SWAP_H_CMC__
#define __COUNT_SWAP_H_CMC__
// Copy from internet https://www.geeksforgeeks.org/number-swaps-sort-adjacent-swapping-allowed/

// C++ program to count number of swaps required 
// to sort an array when only swapping of adjacent 
// elements is allowed. 
#include <bits/stdc++.h> 

/* This function merges two sorted arrays and returns inversion
   count in the arrays.*/
template <typename T>
int merge (vector<T>& arr, vector<T>& temp, int left, int mid, int right)
{
    int inv_count = 0;

    int i = left; /* i is index for left subarray*/
    int j = mid;  /* i is index for right subarray*/
    int k = left; /* i is index for resultant merged subarray*/
    while ((i <= mid - 1) && (j <= right))
    {
        if (arr[i] <= arr[j])
            temp[k++] = arr[i++];
        else
        {
            temp[k++] = arr[j++];
            /* this is tricky -- see above explanation/
              diagram for merge()*/
            inv_count = inv_count + (mid - i);
        }
    }

    /* Copy the remaining elements of left subarray
     (if there are any) to temp*/
    while (i <= mid - 1)
        temp[k++] = arr[i++];

    /* Copy the remaining elements of right subarray
     (if there are any) to temp*/
    while (j <= right)
        temp[k++] = arr[j++];

    /*Copy back the merged elements to original array*/
    for (i=left; i <= right; i++)
        arr[i] = temp[i];

    return inv_count;
}

/* An auxiliary recursive function that sorts the input
   array and returns the number of inversions in the
   array. */
template <typename T>
int _mergeSort(vector<T>& arr, vector<T>& temp, int left, int right)
{
    int mid, inv_count = 0;
    if (right > left)
    {
        /* Divide the array into two parts and call
          _mergeSortAndCountInv() for each of the parts */
        mid = (right + left)/2;

        /* Inversion count will be sum of inversions in
           left-part, right-part and number of inversions
           in merging */
        inv_count  = _mergeSort(arr, temp, left, mid);
        inv_count += _mergeSort(arr, temp, mid+1, right);

        /*Merge the two parts*/
        inv_count += merge(arr, temp, left, mid+1, right);
    }
    return inv_count;
}

/* This function sorts the input array and returns the
   number of inversions in the array */
template <typename T>
inline int sort_countSwaps (vector<T>& arr)
{
    int n = arr.size();
    vector<T> temp (n);
    return _mergeSort (arr, temp, 0, n - 1);
}

template <typename T>
inline int countSwaps (vector<T> arr)
{
    return sort_countSwaps (arr);
}
#endif
