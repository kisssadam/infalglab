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

int elementsInArray[] = { 10, 100, 1000, static_cast<int>(std::pow(10.0, 6.0)) };
size_t elementsInArraySize = sizeof(elementsInArray) / sizeof(int);

int maxValues[] = { 2, 4, 10, 1000, static_cast<int>(std::pow(10.0, 8.0)) };
size_t maxValuesSize = sizeof(maxValues) / sizeof(int);

void fillArrayWithRandomNumbers(size_t size, int *array, int maxValue) {
	for (unsigned int i = 0; i < size; ++i) {
		array[i] = std::rand() % maxValue;
	}
}

int *randomArray(size_t size, int maxValue) {
	int *randomNumbers = new int[size];
	
	fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
	
	return randomNumbers;
}

int findMin(size_t size, int* array) {
	unsigned int i = 0;
	unsigned int j = 1;

	while (j < size) {
		if (array[j] < array[i]) {
			i = j;
		}
		j++;
	}
	
	return i;
}

int findMinRecursively(size_t size, int *array, int firstPosition, int lastPosition) {
	if (firstPosition == lastPosition) {
		return firstPosition;
	} else {
		int k = (firstPosition + lastPosition) / 2;

		int leftResult = findMinRecursively(size, array, firstPosition, k);
		int rightResult = findMinRecursively(size, array, k + 1, lastPosition);
		
		return array[leftResult] <= array[rightResult] ? leftResult : rightResult;
	}
}

int findMinRec(size_t size, int *array) {
	return findMinRecursively(size, array, 0, size - 1);
}

double benchmarkFunction(int (*functionToBenchmark)(size_t, int*), size_t size, int *array) {
	auto start = std::chrono::system_clock::now();
	for (unsigned int k = 0; k < MAX_ITERATION; k++) {
		functionToBenchmark(size, array);
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

	double **sequentialResults = allocateResults(elementsInArraySize, maxValuesSize);
	double **recursiveResults = allocateResults(elementsInArraySize, maxValuesSize);

	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		for (unsigned int j = 0; j < maxValuesSize; j++) {
			size_t size = elementsInArray[i];
			int *randomNumbers = randomArray(size, maxValues[i]);

			sequentialResults[i][j] = benchmarkFunction(findMin, size, randomNumbers);
			recursiveResults[i][j] = benchmarkFunction(findMinRec, size, randomNumbers);

			delete[] randomNumbers;
		}
	}

	printInfo();

	printResults("Sequential minimum search algorithm:", sequentialResults);
	std::cout << std::endl;
	printResults("Recursive minimum search algorithm:", recursiveResults);

	freeResults(sequentialResults);
	freeResults(recursiveResults);

	return EXIT_SUCCESS;
}
