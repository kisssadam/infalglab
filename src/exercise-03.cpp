#include <cstdlib>
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
	for (unsigned int i = 0; i < size; ++i) {
		array[i] = randomInt(maxValue);
	}
}

int *randomArray(size_t size, int maxValue) {
	int *randomNumbers = new int[size];
	
	fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
	
	return randomNumbers;
}

void minimumSelectionOrder(size_t size, int* array) {
	for (unsigned i = 0; i < size - 1; i++) {
		for (unsigned int j = i + 1; j < size; j++) {
			if (array[j] < array[i]) {
				std::swap(array[j], array[i]);
			}
		}
	}
}

int binarySearch(size_t size, int* array, int numberToFind) {
	int lower = 0;
	int upper = size - 1;
	
	while (upper > lower) {
		/*
			További információ arról, hogy miért nem a ((lower + upper) / 2) kifejezést használtam:
			https://en.wikipedia.org/wiki/Binary_search_algorithm
			https://helloacm.com/c-coding-exercise-first-bad-version/
		*/
		int k = lower + (upper - lower) / 2;
		
		if (array[k] == numberToFind) {
			lower = upper = k;
		} else if (array[k] > numberToFind) {
			upper = k - 1;
		} else {
			lower = k + 1;
		}
	}

	return lower <= upper && array[lower] == numberToFind ? lower : size;
}

double benchmarkOrderingAlgorithm(void (*functionToBenchmark)(size_t, int*), size_t size, int maxValue) {
	int *randomNumbers = randomArray(size, maxValue);
	double elapsedSeconds = 0.0;

	for (unsigned int k = 0; k < MAX_ITERATION; k++) {
		fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
		
		auto start = std::chrono::system_clock::now();
		functionToBenchmark(size, randomNumbers);
		auto end = std::chrono::system_clock::now();

		elapsedSeconds += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	}
	
	delete[] randomNumbers;

	return elapsedSeconds / static_cast<double>(MAX_ITERATION);
}

double benchmarkSearchFunction(int (*functionToBenchmark)(size_t, int*, int), size_t size, int *array, int maxValue) {
	auto start = std::chrono::system_clock::now();
	for (unsigned int k = 0; k < MAX_ITERATION; k++) {
		int numberToFind = randomInt(maxValue);
		functionToBenchmark(size, array, numberToFind);
	}
	auto end = std::chrono::system_clock::now();

	double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	return elapsedSeconds / static_cast<double>(MAX_ITERATION);
}

double **allocateResults(size_t rowNum, size_t columnNum) {
	double **results = new double*[rowNum];

	for (unsigned int i = 0; i < rowNum; i++) {
		results[i] = new double[columnNum];
	}

	return results;
}

void freeResults(double **results) {
	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		delete[] results[i];
	}

	delete[] results;
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
	for (unsigned int i = 0; i < maxValuesSize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << maxValues[i];
		std::cout << (i + 1 == maxValuesSize ? "\n" : " ");
	}

	std::cout.precision(PRECISION);
	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		std::cout << std::setw(COLUMN_LENGTH) << std::fixed << elementsInArray[i] << " ";
		
		for (unsigned int j = 0; j < maxValuesSize; j++) {
			std::cout << std::setw(COLUMN_LENGTH) << std::fixed << results[i][j];
			std::cout << (j + 1 == maxValuesSize ? "\n" : " ");
		}
	}
}

int main() {
	std::srand(time(NULL));

	double **orderingResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **binarySearchResults = allocateResults(elementsInArraySize, maxValuesSize);

	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		for (unsigned int j = 0; j < maxValuesSize; j++) {
			size_t size = elementsInArray[i];

			orderingResults[i][j] = benchmarkOrderingAlgorithm(minimumSelectionOrder, size, maxValues[j]);
			
			int *randomNumbers = randomArray(size, maxValues[i]);
			minimumSelectionOrder(size, randomNumbers);
			binarySearchResults[i][j] = benchmarkSearchFunction(binarySearch, size, randomNumbers, maxValues[j]);

			delete[] randomNumbers;
		}
	}

	printInfo();

	printResults("Ordering algorithm:", orderingResults);
	std::cout << std::endl;
	printResults("Binary search algorithm:", binarySearchResults);

	freeResults(orderingResults);
	freeResults(binarySearchResults);

	return EXIT_SUCCESS;
}
