#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>

#define DEVELOPER_NAME "Kiss Sándor Ádám"

#define MAX_ITERATION 10
#define PRECISION 10
#define COLUMN_LENGTH 15

const int sizesToValidateSortFunction[] = {1, 2, 5, 10, 100, 1000};
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

/////////////////////////////////////////////////////////////////////////////////////////////

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

void quickSort(int* array, int lower, int upper) {
	int pivot = array[upper];
	int i = lower;
	int j = upper;

	while (i <= j) {
		while (array[i] < pivot) {
			i++;
		}
		while (array[j] > pivot) {
			j--;
		}
		if (i <= j) {
			std::swap(array[i++], array[j--]);
		}
	}
	
	if (lower < j) {
		quickSort(array, lower, j);
	}
	if (i < upper) {
		quickSort(array, i, upper);
	}
}

void quickSort(int size,  int* array) {
	quickSort(array, 0, size - 1);
}

// leszámláló rendezés
void countingSort(int size, int *array, int maxValue) {
	int* C = new int[maxValue]();
	
	for (int i = 0; i < size; i++) {
		C[array[i]]++;
	}
	
	for (int i = 0, j = 0; i < size; ) {
		if (C[j] > 0) {
			array[i++] = j;
			C[j]--;
		} else {
			j++;
		}
	}

	delete[] C;
}

int bucketKeyFunction(int element) {
	return element;
}

// edényrendezés
void bucketSort(int size, int *array, int maxValue) {
	std::multimap<int, int> *buckets = new std::multimap<int, int>();

	for (int i = 0; i < size; ++i) {
		buckets->emplace(bucketKeyFunction(array[i]), array[i]);
	}

	int j = 0;
	for (std::multimap<int,int>::iterator it = buckets->begin(); it != buckets->end(); ++it) {
		array[j++] = it->second;
	}

	delete buckets;
}

// edényrendezés mutatókkal
void bucketSortWithPointers(int size, int *array, int maxValue) {
	int *buckets = new int[size];
	int **bucketPointers = new int*[maxValue];
	
	for (int i = 0; i < maxValue; ++i) {
		bucketPointers[i] = buckets + size;
	}

	for (int i = size - 1; i >= 0; --i) {
		for (int j = array[i]; j >= 0; --j) {
			--bucketPointers[j];
		}
	}

	for (int i = 0; i < size; ++i) {
		*bucketPointers[array[i]]++ = array[i];
	}

	std::memcpy(array, buckets, size * sizeof(int));

	delete[] buckets;
	delete[] bucketPointers;
}

int bucketHashFunction(int key, int size) {
	return key / size;
}

void bucketSortWithHashFunction(int size, int *array, int maxValue) {
	std::vector<int>** buckets = new std::vector<int>*[maxValue];
	for (int i = 0; i < maxValue; ++i) {
		buckets[i] = new std::vector<int>();
	}

	for (int i = 0; i < size; ++i) {
		buckets[bucketHashFunction(bucketKeyFunction(array[i]), size)]->push_back(array[i]);
	}

	for (int i = 0, j = 0; i < maxValue; ++i) {
		std::sort(buckets[i]->begin(), buckets[i]->end());
		for (std::vector<int>::iterator it = buckets[i]->begin(); it != buckets[i]->end(); ++it) {
			array[j++] = *it;
		}
	}

	for (int i = 0; i < maxValue; ++i) {
		delete buckets[i];
	}
	delete[] buckets;
}

void radixSort(int size, int *array, int maxValue, int g) {
	std::multimap<int, int> *bucketsC = new std::multimap<int, int>();
	std::multimap<int, int> *bucketsD = new std::multimap<int, int>();
	
	int *S = new int[size]();
	std::memcpy(S, array, size * sizeof(int));

	int d = maxValue - 1;
	while (d > 0) {
		d = floor(d / g);

		bucketsC->clear();
		bucketsD->clear();
		
		for (int i = 0; i < size; ++i) {
			bucketsC->emplace(S[i] % g, array[i]);
			bucketsD->emplace(S[i] % g, floor(S[i] / g));
		}

		int i = 0;
		for (std::multimap<int,int>::iterator it_C = bucketsC->begin(), it_D = bucketsD->begin();
				it_C != bucketsC->end() && it_D != bucketsD->end();
				++it_C, ++it_D) {
			array[i] = it_C->second;
			S[i++] = it_D->second;
		}
	}

	delete bucketsC;
	delete bucketsD;
	delete[] S;
}

/////////////////////////////////////////////////////////////////////////////////////////////

double measureExecutionTime(void (*functionToBenchmark)(int, int*, int, int), int size, int *array, int maxValue, int g) {
	auto start = std::chrono::system_clock::now();
	functionToBenchmark(size, array, maxValue, g);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

double benchmarkSortingAlgorithm(void (*functionToBenchmark)(int, int*, int, int), int size, int maxValue, int maxIteration, int g) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);
	
	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, randomNumbers, maxValue, g);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

/////////////////////////////////////////////////////////////////////////////////////////////

double measureExecutionTime(void (*functionToBenchmark)(int, int*, int), int size, int *array, int maxValue) {
	auto start = std::chrono::system_clock::now();
	functionToBenchmark(size, array, maxValue);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

double benchmarkSortingAlgorithm(void (*functionToBenchmark)(int, int*, int), int size, int maxValue, int maxIteration) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);
	
	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, randomNumbers, maxValue);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

/////////////////////////////////////////////////////////////////////////////////////////////

double measureExecutionTime(int* (*functionToBenchmark)(int, int*), int size, int *array) {
	auto start = std::chrono::system_clock::now();
	delete[] functionToBenchmark(size, array);
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

/////////////////////////////////////////////////////////////////////////////////////////////

double measureExecutionTime(void (*functionToBenchmark)(int, int*), int size, int *array) {
	auto start = std::chrono::system_clock::now();
	functionToBenchmark(size, array);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
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

/////////////////////////////////////////////////////////////////////////////////////////////

double benchmarkHeapSortingAlgorithm(void (*functionToBenchmark)(int, int*),
									 void (*heapifyFunction)(int, int*),
									 int size,
									 int maxValue,
									 int maxIteration) {
	double elapsedSeconds = 0.0;
	
	for (int k = 0; k < maxIteration; k++) {
		int *heap = randomHeap(heapifyFunction, size, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, heap);
		delete[] heap;
	}

	return elapsedSeconds / static_cast<double>(maxIteration);
}

/////////////////////////////////////////////////////////////////////////////////////////////

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

bool isSortFunctionValid(void (*sortFunctionToValidate)(int, int*, int)) {
	bool isFunctionValid = true;

	for (int i = 0; i < lengthOfSizesToValidateSortFunction; i++) {
		int size = sizesToValidateSortFunction[i];
		int *descendingNumbers = descendingNumberArray(size);
		
		const int maxValue = size + 1;
		sortFunctionToValidate(size, descendingNumbers, maxValue);
		
		if (isArraySorted(size, descendingNumbers) == false) {
			isFunctionValid = false;
			break;
		}
		
		delete[] descendingNumbers;
	}
	
	return isFunctionValid;
}

bool isSortFunctionValid(void (*sortFunctionToValidate)(int, int*, int, int)) {
	bool isFunctionValid = true;

	for (int i = 0; i < lengthOfSizesToValidateSortFunction; ++i) {
		for (int g = 2; g <= 20; ++g) {
			int size = sizesToValidateSortFunction[i];
			int *descendingNumbers = descendingNumberArray(size);
			
			const int maxValue = size + 1;
			sortFunctionToValidate(size, descendingNumbers, maxValue, g);
			
			if (isArraySorted(size, descendingNumbers) == false) {
				isFunctionValid = false;
				break;
			}
			
			delete[] descendingNumbers;
		}
	}
	
	return isFunctionValid;
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

/////////////////////////////////////////////////////////////////////////////////////////////

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
	
	std::cout.precision(PRECISION);

	{
		const int size = 100000;
		const int maxValue = 100000000;
		
		double elapsedSeconds = benchmarkSortingAlgorithm(countingSort, size, maxValue, MAX_ITERATION);
		std::cout << std::setw(10) << std::fixed << "Leszámláló rendezés\telapsedSeconds: " << elapsedSeconds << "\tsize: " << size << "\tmaxValue: " << maxValue << std::endl;
	}
	std::cout << std::endl;

	for (int size = 1; size <= 10000; size *= 10) {
		int maxValue = size * size;
		double elapsedSeconds = benchmarkSortingAlgorithm(bucketSort, size, maxValue, MAX_ITERATION);
		std::cout << std::setw(10) << std::fixed << "Edényrendezés\telapsedSeconds: " << elapsedSeconds << "\tsize: " << size << "\tmaxValue: " << maxValue << std::endl;
	}
	std::cout << std::endl;

	for (int size = 1; size <= 10000; size *= 10) {
		int maxValue = size * size;
		double elapsedSeconds = benchmarkSortingAlgorithm(bucketSortWithHashFunction, size, maxValue, MAX_ITERATION);
		std::cout << std::setw(10) << std::fixed << "Edényrendezés hashtáblával\telapsedSeconds: " << elapsedSeconds << "\tsize: " << size << "\tmaxValue: " << maxValue << std::endl;
	}
	std::cout << std::endl;

	for (int size = 1; size <= 1000; size *= 10) {
		int maxValue = size * size;
		double elapsedSeconds = benchmarkSortingAlgorithm(bucketSortWithPointers, size, maxValue, MAX_ITERATION);
		std::cout << std::setw(10) << std::fixed << "Edényrendezés mutatókkal\telapsedSeconds: " << elapsedSeconds << "\tsize: " << size << "\tmaxValue: " << maxValue << std::endl;
	}
	std::cout << std::endl;

	for (int g = 2; g <= 20; ++g) {
		const int size = 100000;
		const int maxValue = 100000000;
		
		double elapsedSeconds = benchmarkSortingAlgorithm(radixSort, size, maxValue, MAX_ITERATION, g);
		std::cout << std::setw(10) << std::fixed << "Radix rendezés\tg: " << g << "\telapsedSeconds: " << elapsedSeconds << "\tsize: " << size << "\tmaxValue: " << maxValue << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Leszámláló rendezés: " << (isSortFunctionValid(countingSort) == true ? "valid" : "invalid") << std::endl;
	std::cout << "Edényrendezés: " << (isSortFunctionValid(bucketSort) == true ? "valid" : "invalid") << std::endl;
	std::cout << "Edényrendezés mutatókkal: " << (isSortFunctionValid(bucketSortWithPointers) == true ? "valid" : "invalid") << std::endl;
	std::cout << "Edényrendezés hashtáblával: " << (isSortFunctionValid(bucketSortWithHashFunction) == true ? "valid" : "invalid") << std::endl;
	std::cout << "Radix rendezés: " << (isSortFunctionValid(radixSort) == true ? "valid" : "invalid") << std::endl;

	return EXIT_SUCCESS;
}
