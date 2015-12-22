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

const int elementsInArray[] = { 1000, 10000, 100000 };
const size_t elementsInArraySize = sizeof(elementsInArray) / sizeof(int);

const int maxValues[] = { 2, 4, 10, 1000 };
const size_t maxValuesSize = sizeof(maxValues) / sizeof(int);

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

void insertionSort(size_t size, int* array) {
	for (size_t i = 1; i < size; i++) {
		int temp = array[i];
		size_t j = i;
		
		while ((j > 0) and (temp < array[j - 1])) {
			array[j] = array[j - 1];
			j--;
		}
		
		array[j] = temp;
	}
}

void shellSort(size_t size, int* array) {
	const int k[] = {701, 301, 132, 57, 23, 10, 4, 1};
	const int l = sizeof(k) / sizeof(int);
	
	for (int s =0; s < l; s++) {
		for (size_t i = k[s]; i < size; i++) {
			int j = i;
			int temp = array[j];
			
			while ((j >= k[s]) and (array[j - k[s]] > temp)) {
				array[j] = array[j - k[s]];
				j = j - k[s];
			}
			
			array[j] = temp;
		}
	}
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

	if (size >= 2) {
		for (size_t i = 0; i < (size - 2); i++) {
			if (array[i] > array[i + 1]) {
				isSorted = false;
				break;
			}
		}
	}

	return isSorted;
}

bool isSortFunctionValid(void (*sortFunctionToValidate)(size_t, int*)) {
	const int sizes[] = {1, 2, 5, 10, 100, 1000, 10000};
	const int lengthOfSizes = sizeof(sizes) / sizeof(int);

	bool isFunctionValid = true;

	for (size_t i = 0; i < lengthOfSizes; i++) {
		size_t size = sizes[i];
		
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

	double **insertionSortResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **shellSortResults = allocateResults(elementsInArraySize, maxValuesSize);

	for (size_t i = 0; i < elementsInArraySize; i++) {
		for (size_t j = 0; j < maxValuesSize; j++) {
			size_t size = elementsInArray[i];
			int maxValue = maxValues[j];

			insertionSortResults[i][j] = benchmarkSortingAlgorithm(insertionSort, size, maxValue, MAX_ITERATION);
			shellSortResults[i][j] = benchmarkSortingAlgorithm(shellSort, size, maxValue, MAX_ITERATION);
		}
	}

	printInfo();
	printResults("Insertion sort algorithm:", insertionSortResults);
	printResults("Shell sort algorithm:", shellSortResults);

	freeResults(insertionSortResults);
	freeResults(shellSortResults);

	bool isInsertionSortFunctionValid = isSortFunctionValid(insertionSort);
	bool isShellSortFunctionValid = isSortFunctionValid(shellSort);
	
	std::cout << "Insertion sort function: " << (isInsertionSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << "Shell sort function: " << (isShellSortFunctionValid == true ? "valid" : "invalid") << std::endl;
	std::cout << std::endl;

	return EXIT_SUCCESS;
}
