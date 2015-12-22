# infalglab
Programs created for the "Computer Algorithms" subject at University of Debrecen.
These programs contain many implemented search and sort algorithm. Every program will print out the average execution time of each implemented algorithm.

## Programs with algorithms
* exercise-01: iterative search algorithm
* exercise-02: iterative and recursive minimum search algorithm
* exercise-02-multithreaded: multithreaded recursive minimum search algorithm in Java
* exercise-03: minimum selection sort, binary search
* exercise-04: bubble sort, merge sort
* exercise-05: insertion sort, shell sort
* exercise-06: matrix sort
* exercise-07: heapsort
* exercise-08: heapsort with linear heapify function
* exercise-09: quicksort (pivot is the last element of the array, but it should be: (array[lower] + array[upper]) / 2))
* exercise-10: counting sort, bucket sort, bucket sort with hashtable, radix sort
* exercise-11: binary search tree
* exercise-written-examination: I had to write a function that takes two parameters: an array and the size of the array. The function had to return two pointers to two different arrays. The elements in the left array had to be lower or equal to the elements in the other one.

## How to compile them?
Before compiling the c++ programs you will need GNU Make and g++ with c++11 support. At this time my compiler version is ```g++ (GCC) 5.3.1 20151207 (Red Hat 5.3.1-2)```. To compile them, please type ```make``` into your terminal.

To compile the exercise-02-multithreaded program you will need a Java 8 compiler. Just simply type ```javac MinimumSearch.java``` into your terminal.
