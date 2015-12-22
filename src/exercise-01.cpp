#include <cstdlib>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

#define DEVELOPER_NAME "Kiss Sándor Ádám"

#define MAX_ITERATION 100
#define PRECISION 10
#define COLUMN_LENGTH 15

int elementsInArray[] = { 10, 100, 1000, 10000, static_cast<int>(std::pow(10.0, 6.0)), static_cast<int>(std::pow(10.0, 8.0)) };
size_t elementsInArraySize = sizeof(elementsInArray) / sizeof(int);

int maxValues[] = { 2, 4, 10, 1000, 10000, static_cast<int>(std::pow(10.0, 8.0)), static_cast<int>(std::pow(10.0, 9.0)) };
size_t maxValuesSize = sizeof(maxValues) / sizeof(int);

double **results;

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

int find(size_t size, int *array, int valueToFind) {
	unsigned int i = 0;

	while (i < size && array[i] != valueToFind) {
		i++;
	}

	return i;
}

double **allocateResults(size_t rowNum, size_t columnNum) {
	double **results = new double*[rowNum];

	for (unsigned int i = 0; i < rowNum; i++) {
		results[i] = new double[columnNum];
	}

	return results;
}

void freeResults() {
	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		delete[] results[i];
	}

	delete[] results;
}

void printResults() {
	std::cout << "Developed by " << DEVELOPER_NAME << std::endl;
	std::cout << std::endl;
	std::cout << "n: the number of elements in the generated array" << std::endl;
	std::cout << "m: the maximum possible value in the array" << std::endl;
	std::cout << "Results are in seconds!" << std::endl;
	std::cout << std::endl;

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

	results = allocateResults(elementsInArraySize, maxValuesSize);

	for (unsigned int i = 0; i < elementsInArraySize; i++) {
		for (unsigned int j = 0; j < maxValuesSize; j++) {
			size_t size = elementsInArray[i];
			
			int *randomNumbers = randomArray(size, maxValues[i]);

			auto start = std::chrono::system_clock::now();
			for (unsigned int k = 0; k < MAX_ITERATION; k++) {
				int valueToFind = std::rand() % maxValues[i];
				find(size, randomNumbers, valueToFind);
			}
			auto end = std::chrono::system_clock::now();

			/*
				// good for milliseconds:
				auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				std::cout << (double) elapsed.count() / static_cast<double>(MAX_ITERATION) << std::endl;
			*/

			// seconds
			double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			results[i][j] = elapsed_seconds / static_cast<double>(MAX_ITERATION);

			delete[] randomNumbers;
		}
	}

	printResults();

	freeResults();

	return EXIT_SUCCESS;
}
