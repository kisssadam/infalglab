#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

#define DEVELOPER_NAME "Kiss Sándor Ádám"

#define MAX_ITERATION 1000
#define PRECISION 10
#define COLUMN_LENGTH 15

int elementsInArray[] = { 10, 100, 1000 };
size_t elementsInArraySize = sizeof(elementsInArray) / sizeof(int);

int maxValues[] = { 2, 4, 10, 1000 };
size_t maxValuesSize = sizeof(maxValues) / sizeof(int);

int randomInt(int maxValue) {
	return std::rand() % maxValue;
}

void fillArrayWithRandomNumbers(size_t size, int *array, int maxValue) {
	for (size_t i = 0; i < size; ++i) {
		array[i] = randomInt(maxValue);
	}
}

int *randomArray(size_t size, int maxValue) {
	int *randomNumbers = new int[size];
	
	fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
	
	return randomNumbers;
}

void fillArrayWithDescendingNumbers(size_t size, int *array) {
	for (size_t i = 0, value = size; i < size; i++, value--) {
		array[i] = value;
	}
}

int *descendingNumberArray(size_t size) {
	int *descendingNumbers = new int[size];

	fillArrayWithDescendingNumbers(size, descendingNumbers);

	return descendingNumbers;
}

double **allocateResults(size_t rowNum, size_t columnNum) {
	double **results = new double*[rowNum];

	for (size_t i = 0; i < rowNum; i++) {
		results[i] = new double[columnNum];
	}

	return results;
}

void freeResults(double **results) {
	for (size_t i = 0; i < elementsInArraySize; i++) {
		delete[] results[i];
	}

	delete[] results;
}

void bubbleSort(size_t size, int* array) {
	for (size_t i = 0; i < (size - 2); i++) {
		for (size_t j = (size - 1); i < j; j--) {
			if (array[j] < array[j - 1]) {
				std::swap(array[j], array[j - 1]);
			}
		}
	}
}

void fasterBubbleSort(size_t size, int* array) {
	for (size_t i = 0; i < (size - 2); i++) {
		size_t swapIndex = i;
		
		for (size_t j = (size - 1); i < j; j--) {
			if (array[j] < array[j - 1]) {
				std::swap(array[j], array[j - 1]);
				swapIndex = j - 1;
			}
		}
		
		i = swapIndex;
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

void mergeSort(size_t size, int* array) {
	return mergeSortWithBounds(0, size - 1, array);
}

void hybridMergeSortWithBounds(int lower, int upper, int* array, int x) {
	if (upper-lower <= x) {
		for (int i = 0; i < (upper - 1); i++) {
			int swapIndex = i;
			
			for (int j = (upper); i < j; j--) {
				if (array[j] < array[j - 1]) {
					std::swap(array[j], array[j - 1]);
					swapIndex = j - 1;
				}
			}
			
			i = swapIndex;
		}
	} else {
		int k = lower + (upper - lower) / 2;
		hybridMergeSortWithBounds(lower, k, array, x);
		hybridMergeSortWithBounds(k + 1, upper, array, x);
		
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

void hybridMergeSort(size_t size, int* array) {
	return hybridMergeSortWithBounds(0, size - 1, array, 15);
}

double measureExecutionTimeOfHybridMergeSortingAlgorithm(void (*hybridMergeSortToMeasure)(int, int, int*, int), size_t size, int *array, int x) {
	auto start = std::chrono::system_clock::now();
	hybridMergeSortToMeasure(0, size - 1, array, x);
	auto end = std::chrono::system_clock::now();

	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();	
}

double benchmarkHybridMergeSortingAlgorithm(size_t size, int maxValue, int x, int maxIteration) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);

	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTimeOfHybridMergeSortingAlgorithm(hybridMergeSortWithBounds, size, randomNumbers, x);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

double measureExecutionTime(void (*functionToBenchmark)(size_t, int*), size_t size, int *array) {
	auto start = std::chrono::system_clock::now();
	functionToBenchmark(size, array);
	auto end = std::chrono::system_clock::now();
	
	return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

double benchmarkSortingAlgorithm(void (*functionToBenchmark)(size_t, int*), size_t size, int maxValue, int maxIteration) {
	double elapsedSeconds = 0.0;
	int *randomNumbers = randomArray(size, maxValue);
	
	for (int k = 0; k < maxIteration; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		elapsedSeconds += measureExecutionTime(functionToBenchmark, size, randomNumbers);
	}

	delete[] randomNumbers;
	return elapsedSeconds / static_cast<double>(maxIteration);
}

bool isArraySorted(size_t size, int *array) {
	bool isSorted = true;

	for (size_t i = 0; i < (size - 2); i++) {
		if (array[i] > array[i + 1]) {
			isSorted = false;
			break;
		}
	}

	return isSorted;
}

bool isSortFunctionValid(void (*sortFunctionToValidate)(size_t, int*)) {
	bool isFunctionValid = true;

	for (size_t i = 0; i < elementsInArraySize; i++) {
		size_t size = elementsInArray[i];
		
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

void determineOptimalXForHybridMergeSort() {
	const size_t size = 10000;
	const int maxValue = 10000;

	std::cout << "Calculating the optimal x for hybrid merge sort. Array size fixed: " << size << std::endl;

	double fasterBubbleSortResult = benchmarkSortingAlgorithm(fasterBubbleSort, size, maxValue, 10);
	
	std::cout.precision(PRECISION);
	std::cout << std::fixed << "faster bubble sort: " << fasterBubbleSortResult << std::endl;

	int startedToRunFasterIndex = -1;
	for (int x = 0; (startedToRunFasterIndex < 0) or (x < startedToRunFasterIndex + 20); x++) {
		double hybridMergeSortResult = benchmarkHybridMergeSortingAlgorithm(size, maxValue, x, 1);

		std::cout << std::fixed << "x: " << x << "\thybrid merge sort: " << hybridMergeSortResult;
		if (hybridMergeSortResult < fasterBubbleSortResult) {
			if (startedToRunFasterIndex < 0) {
				startedToRunFasterIndex = x;
			}
			std::cout << " <--------";
		}
		std::cout << std::endl;
	}

	std::cout << "Merge sort is faster than bubble sort using x: " << startedToRunFasterIndex << std::endl;
	std::cout << std::endl;
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
	for (size_t i = 0; i < maxValuesSize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << maxValues[i];
		std::cout << (i + 1 == maxValuesSize ? "\n" : " ");
	}

	std::cout.precision(PRECISION);
	for (size_t i = 0; i < elementsInArraySize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << elementsInArray[i] << " ";
		
		for (size_t j = 0; j < maxValuesSize; j++) {
			std::cout << std::setw(COLUMN_LENGTH) << std::fixed << results[i][j];
			std::cout << (j + 1 == maxValuesSize ? "\n" : " ");
		}
	}
	std::cout << std::endl;
}

int main(void) {
	std::srand(time(NULL));

	double **bubbleSortResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **fasterBubbleSortResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **mergeSortResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **hybridMergeSortResults = allocateResults(elementsInArraySize, maxValuesSize);

	for (size_t i = 0; i < elementsInArraySize; i++) {
		for (size_t j = 0; j < maxValuesSize; j++) {
			size_t size = elementsInArray[i];
			int maxValue = maxValues[j];

			bubbleSortResults[i][j] = benchmarkSortingAlgorithm(bubbleSort, size, maxValue, MAX_ITERATION);
			fasterBubbleSortResults[i][j] = benchmarkSortingAlgorithm(fasterBubbleSort, size, maxValue, MAX_ITERATION);
			mergeSortResults[i][j] = benchmarkSortingAlgorithm(mergeSort, size, maxValue, MAX_ITERATION);
			hybridMergeSortResults[i][j] = benchmarkSortingAlgorithm(hybridMergeSort, size, maxValue, MAX_ITERATION);
		}
	}

	printInfo();
	printResults("Bubble sort algorithm:", bubbleSortResults);
	printResults("Faster bubble sort algorithm:", fasterBubbleSortResults);
	printResults("Merge sort algorithm:", mergeSortResults);
	printResults("Hybrid merge sort algorithm:", hybridMergeSortResults);

	freeResults(bubbleSortResults);
	freeResults(fasterBubbleSortResults);
	freeResults(mergeSortResults);
	freeResults(hybridMergeSortResults);

	bool isBubbleSortFunctionValid = isSortFunctionValid(bubbleSort);
	bool isFasterBubbleSortFunctionValid = isSortFunctionValid(fasterBubbleSort);
	bool isMergeSortFunctionValid = isSortFunctionValid(mergeSort);
	bool isHybridMergeSortFunctionValid = isSortFunctionValid(mergeSort);

	std::cout << "Bubble sort function: " << (isBubbleSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << "Faster bubble sort function: " << (isFasterBubbleSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << "Merge sort function: " << (isMergeSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << "Hybrid merge sort function: " << (isHybridMergeSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << std::endl;

	determineOptimalXForHybridMergeSort();

	return EXIT_SUCCESS;
}
