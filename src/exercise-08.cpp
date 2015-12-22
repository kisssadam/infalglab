#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

#define DEVELOPER_NAME "Kiss Sándor Ádám"

#define MAX_ITERATION 10
#define PRECISION 10
#define COLUMN_LENGTH 15


const int sizesToValidateSortFunction[] = {1, 2, 5, 10, 100, 1000, 10000};
const int lengthOfSizesToValidateSortFunction = sizeof(sizesToValidateSortFunction) / sizeof(int);

const int elementsInArray[] = { 10, 100, 1000, 10000, 100000, 1000000 };
const int elementsInArraySize = sizeof(elementsInArray) / sizeof(int);

const int maxValues[] = { 2, 4, 10, 1000 };
const int maxValuesSize = sizeof(maxValues) / sizeof(int);

int randomInt(int maxValue) {
	return std::rand() % maxValue;
}

void fillArrayWithRandomNumbers(int size, int *array, int maxValue) {
	for (int i = 0; i < size; ++i) {
		array[i] = randomInt(maxValue);
	}
}

int *randomArray(int size, int maxValue) {
	int *randomNumbers = new int[size];
	
	fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
	
	return randomNumbers;
}

void fillArrayWithDescendingNumbers(int size, int *array) {
	for (int i = 0, value = size; i < size; i++, value--) {
		array[i] = value;
	}
}

int *descendingNumberArray(int size) {
	int *descendingNumbers = new int[size];

	fillArrayWithDescendingNumbers(size, descendingNumbers);

	return descendingNumbers;
}

double **allocateResults(int rowNum, int columnNum) {
	double **results = new double*[rowNum];

	for (int i = 0; i < rowNum; i++) {
		results[i] = new double[columnNum];
	}

	return results;
}

void freeResults(double **results) {
	for (int i = 0; i < elementsInArraySize; i++) {
		delete[] results[i];
	}

	delete[] results;
}

void insertionSort(int size, int* array) {
	for (int i = 1; i < size; i++) {
		int temp = array[i];
		int j = i;
		
		while ((j > 0) and (temp < array[j - 1])) {
			array[j] = array[j - 1];
			j--;
		}
		
		array[j] = temp;
	}
}

int* matrixSort(int size, int* array) {
	int *B = new int[size]();
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < i; j++) {
			B[array[j] <= array[i] ? i : j]++;
		}
	}
	
	int* result = new int[size];
	for (int i = 0; i < size; i++) {
		result[B[i]] = array[i];
	}
	delete[] B;

	return result;
}

void heapify(int size, int* array) {
	for (int i = 1; i <= size - 1; i++) {
		int j = i;
		while (j > 0 and array[(i - 1) / 2] < array[j]) {
			std::swap(array[j], array[(i - 1) / 2]);
			j = (i - 1) / 2;
		}
	}
}

void heapifyLinear(int size, int* array) {
	for (int i = (size / 2) - 1; i > 0; i--) {
		for (int j = i; j <= (size / 2) - 1; j--) {
			int k = 2 * j + 1;
			if (k < size - 1 and array[k + 1] > array[k]) {
				k++;
			}
			if (array[j] < array[k]) {
				std::swap(array[j], array[k]);
				j = k;
			} else {
				j = size;
			}
		}
	}
}

int* randomHeap(void (*heapifyFunction)(int, int*), int size, int maxValue) {
	int* array = randomArray(size, maxValue);
	heapifyFunction(size, array);
	return array;
}

void heapSort(int size, int* array) {
	for (int i = size - 1; i >= 1; i--) {
		std::swap(array[i], array[0]);
		int j = 0;
		while (2 * j + 1 <= i - 1) {
			int k = 2 * j + 1;
			if (k < i - 1 and array[k] < array[k+1]) {
				k = k + 1;
			}
			if (array[k] > array[j]) {
				std::swap(array[k], array[j]);
				j = k;
			} else {
				j = size;
			}
		}
	}
}

void mergeSortWithBounds(int lower, int upper, int* array) {
	if (upper-lower <= 1) {
		if (array[upper] < array[lower]) {
			std::swap(array[upper], array[lower]);
		}
	} else {
		int k = lower + (upper - lower) / 2;
		mergeSortWithBounds(lower, k, array);
		mergeSortWithBounds(k + 1, upper, array);
		
		int* sortedArray = new int[upper + 1];
		
		int i = lower;
		int j = k + 1;
		
		while ((i <= k) and (j <= upper)) {
			int index = i + j - (lower + k + 1);
			sortedArray[index] = (array[i] < array[j] ? array[i++] : array[j++]);
		}
		
		if (i <= k) {
			std::memmove(&array[upper-(k-i)], &array[i], (k-i+1)*sizeof(int));
		}
		std::memcpy(&array[lower], sortedArray, (i+j-(lower+k+1)-1+1)*sizeof(int));

		delete[] sortedArray;
	}
}

void mergeSort(int size, int* array) {
	return mergeSortWithBounds(0, size - 1, array);
}

double measureExecutionTime(int* (*functionToBenchmark)(int, int*), int size, int *array) {
	auto start = std::chrono::system_clock::now();
	delete[] functionToBenchmark(size, array);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

double measureExecutionTime(void (*functionToBenchmark)(int, int*), int size, int *array) {
	auto start = std::chrono::system_clock::now();
	functionToBenchmark(size, array);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

double benchmarkSortingAlgorithm(int* (*functionToBenchmark)(int, int*), int size, int maxValue, int maxIteration) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);
	
	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, randomNumbers);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

double benchmarkSortingAlgorithm(void (*functionToBenchmark)(int, int*), int size, int maxValue, int maxIteration) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);
	
	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, randomNumbers);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

double benchmarkHeapSortingAlgorithm(void (*functionToBenchmark)(int, int*), void (*heapifyFunction)(int, int*), int size, int maxValue, int maxIteration) {
	double elapsedSeconds = 0.0;
	
	for (int k = 0; k < maxIteration; k++) {
		int *heap = randomHeap(heapifyFunction, size, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, heap);
		delete[] heap;
	}

	return elapsedSeconds / static_cast<double>(maxIteration);
}

bool isArraySorted(int size, int *array) {
	bool isSorted = true;

	if (size >= 2) {
		for (int i = 0; i < (size - 2); i++) {
			if (array[i] > array[i + 1]) {
				isSorted = false;
				break;
			}
		}
	}

	return isSorted;
}

bool isSortFunctionValid(void (*sortFunctionToValidate)(int, int*)) {
	bool isFunctionValid = true;

	for (int i = 0; i < lengthOfSizesToValidateSortFunction; i++) {
		int size = sizesToValidateSortFunction[i];
		int *descendingNumbers = descendingNumberArray(size);
		
		sortFunctionToValidate(size, descendingNumbers);
		
		if (isArraySorted(size, descendingNumbers) == false) {
			isFunctionValid = false;
			break;
		}
		
		delete[] descendingNumbers;
	}
	
	return isFunctionValid;
}

bool isSortFunctionValid(int *(*sortFunctionToValidate)(int, int*)) {
	bool isFunctionValid = true;

	for (int i = 0; i < lengthOfSizesToValidateSortFunction; i++) {
		int size = sizesToValidateSortFunction[i];
		int *descendingNumbers = descendingNumberArray(size);
		
		int *result = sortFunctionToValidate(size, descendingNumbers);
		
		if (isArraySorted(size, result) == false) {
			isFunctionValid = false;
			break;
		}
		
		delete[] result;
		delete[] descendingNumbers;
	}
	
	return isFunctionValid;
}

void printInfo(void) {
	std::cout << "Developed by " << DEVELOPER_NAME << std::endl;
	std::cout << std::endl;
	std::cout << "n: the number of elements in the generated array" << std::endl;
	std::cout << "m: the maximum possible value in the array" << std::endl;
	std::cout << "Results are in seconds!" << std::endl;
	std::cout << std::endl;
}

void printResults(std::string title, double **results) {
	std::cout << title << std::endl;
	std::cout << std::setw(COLUMN_LENGTH) << "n \\ m";
	std::cout << " ";
	for (int i = 0; i < maxValuesSize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << maxValues[i];
		std::cout << (i + 1 == maxValuesSize ? "\n" : " ");
	}

	std::cout.precision(PRECISION);
	for (int i = 0; i < elementsInArraySize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << elementsInArray[i] << " ";
		
		for (int j = 0; j < maxValuesSize; j++) {
			std::cout << std::setw(COLUMN_LENGTH) << std::fixed << results[i][j];
			std::cout << (j + 1 == maxValuesSize ? "\n" : " ");
		}
	}
	std::cout << std::endl;
}

int main(void) {
	std::srand(time(NULL));

	double **heapSortResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **heapSortLinearHeapifyResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **mergeSortResults = allocateResults(elementsInArraySize, maxValuesSize);

	for (int i = 0; i < elementsInArraySize; i++) {
		for (int j = 0; j < maxValuesSize; j++) {
			int size = elementsInArray[i];
			int maxValue = maxValues[j];

			heapSortResults[i][j] = benchmarkHeapSortingAlgorithm(heapSort, heapify, size, maxValue, MAX_ITERATION);
			heapSortLinearHeapifyResults[i][j] = benchmarkHeapSortingAlgorithm(heapSort, heapifyLinear, size, maxValue, MAX_ITERATION);
			mergeSortResults[i][j] = benchmarkSortingAlgorithm(mergeSort, size, maxValue, MAX_ITERATION);
		}
	}

	printInfo();
	printResults("Heap sort algorithm:", heapSortResults);
	printResults("Heap sort with linear heapify function algorithm:", heapSortLinearHeapifyResults);
	printResults("Merge sort algorithm:", mergeSortResults);

	freeResults(heapSortResults);
	freeResults(heapSortLinearHeapifyResults);
	freeResults(mergeSortResults);

	return EXIT_SUCCESS;
}
