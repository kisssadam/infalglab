import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.function.Function;

/**
 * 
 * @author Kiss Sándor Ádám
 *
 */
public class MinimumSearch {

	private static final int[] elementsInArray = { 10, 100, 1000, (int) Math.pow(10, 6) };
	private static final int[] maxValues = { 2, 4, 10, 1000, (int) Math.pow(10, 8) };
	private static final Random random = new Random();
	private static final int MAX_ITERATION = 100;

	private static final ForkJoinPool forkJoinPool = new ForkJoinPool();
	private static final Function<int[], Integer> sequentialMinimumSearch = numbers -> sequentialMinimumSearch(numbers);
	private static final Function<int[], Integer> recursiveMinimumSearch = numbers -> recursiveMinimumSearch(numbers);
	private static final Function<int[], Integer> parallelMinimumSearch = numbers -> parallelMinimumSearch(numbers);

	private static final String[] HEADER = { "Size of array", "Max value", "Sequential", "Recursive", "Multithreaded" };
	private static final String CORRECT_MESSAGE = "The program seems to be correct!";
	private static final String NOT_CORRECT_MESSAGE = "It looks like the program is not correct!";

	private static final int LEADING_SPACES = calculateLeadingSpaces();
	private static final int COLUMN_WIDTH = 15;

	private static final int calculateLeadingSpaces() {
		int requiredDigits = (int) (Math.log10(elementsInArray.length * maxValues.length) + 1);
		return requiredDigits + 2;
	}

	public static void main(String[] args) {
		printHeader(LEADING_SPACES, COLUMN_WIDTH, HEADER);

		int outputRowCounter = 0;
		for (int elementInArray : elementsInArray) {
			for (int maxValue : maxValues) {
				printFirstRowPart(++outputRowCounter);

				int[] randomNumbers = random.ints(elementInArray, 0, maxValue).toArray();

				long sequentialExecutionTime = calculateExecutionTime(randomNumbers, sequentialMinimumSearch);
				long recursiveExecutionTime = calculateExecutionTime(randomNumbers, recursiveMinimumSearch);
				long multithreadedExecutionTime = calculateExecutionTime(randomNumbers, parallelMinimumSearch);

				printLastRowPart(COLUMN_WIDTH, elementInArray, maxValue, sequentialExecutionTime,
						recursiveExecutionTime, multithreadedExecutionTime);
			}
			System.out.println();
		}

		System.out.println("Checking whether the program is correct or not.");
		System.out.println(isProgramCorrect() ? CORRECT_MESSAGE : NOT_CORRECT_MESSAGE);
	}

	private static void printHeader(int leadingSpaces, int columnWidth, String... titles) {
		if (leadingSpaces > 0) {
			System.out.printf("%" + leadingSpaces + "c", ' ');
		}

		for (int i = 0; i < titles.length; i++) {
			String title = titles[i];
			System.out.printf("%" + columnWidth + "s%c", title, (i + 1 == titles.length ? '\n' : '\t'));
		}
	}

	private static void printFirstRowPart(int index) {
		System.out.printf("%" + (LEADING_SPACES - 2) + "d. ", index);
	}

	private static void printLastRowPart(int columnWidth, long... values) {
		for (int i = 0; i < values.length; i++) {
			long value = values[i];
			System.out.printf("%" + columnWidth + "d%c", value, (i + 1 == values.length ? '\n' : '\t'));
		}
	}

	private static long calculateExecutionTime(int[] numbers, Function<int[], Integer> minimumSearchFunction) {
		System.gc();

		long startTime = System.currentTimeMillis();
		for (int i = 0; i < MAX_ITERATION; i++) {
			minimumSearchFunction.apply(numbers);
		}
		long endTime = System.currentTimeMillis();

		double elapsedTime = endTime - startTime;
		return Math.round(elapsedTime / MAX_ITERATION);
	}

	public static int parallelMinimumSearch(int[] numbers) {
		MultithreadedMinimumSearch minimumSearch = new MultithreadedMinimumSearch(numbers);
		return forkJoinPool.invoke(minimumSearch);
	}

	private static class MultithreadedMinimumSearch extends RecursiveTask<Integer> {

		private static final long serialVersionUID = -7794083474094004697L;
		private static final int THRESHOLD = 100;

		private int[] numbers;
		private int inclusiveBeginIndex;
		private int exclusiveEndIndex;

		public MultithreadedMinimumSearch(int[] numbers) {
			this(numbers, 0, numbers.length);
		}

		private MultithreadedMinimumSearch(int[] numbers, int inclusiveBeginIndex, int exclusiveEndIndex) {
			this.numbers = numbers;
			this.inclusiveBeginIndex = inclusiveBeginIndex;
			this.exclusiveEndIndex = exclusiveEndIndex;
		}

		@Override
		protected Integer compute() {
			int length = exclusiveEndIndex - inclusiveBeginIndex;

			if (length < THRESHOLD) {
				return recursiveMinimumSearch(numbers, inclusiveBeginIndex, exclusiveEndIndex - 1);
			}

			int split = length / 2;

			MultithreadedMinimumSearch left = new MultithreadedMinimumSearch(numbers, inclusiveBeginIndex,
					inclusiveBeginIndex + split);
			left.fork();
			MultithreadedMinimumSearch right = new MultithreadedMinimumSearch(numbers, inclusiveBeginIndex + split,
					exclusiveEndIndex);

			int rightIndex = right.compute();
			int leftIndex = left.join();

			return numbers[leftIndex] <= numbers[rightIndex] ? leftIndex : rightIndex;
		}

	}

	public static int recursiveMinimumSearch(int[] numbers) {
		return recursiveMinimumSearch(numbers, 0, numbers.length - 1);
	}

	private static int recursiveMinimumSearch(int[] numbers, int begin, int end) {
		if (begin == end) {
			return begin;
		} else {
			int k = (begin + end) / 2;

			int leftIndex = recursiveMinimumSearch(numbers, begin, k);
			int rightIndex = recursiveMinimumSearch(numbers, k + 1, end);

			return numbers[leftIndex] <= numbers[rightIndex] ? leftIndex : rightIndex;
		}
	}

	public static int sequentialMinimumSearch(int[] numbers) {
		int i = 0;

		for (int j = 1; j < numbers.length; j++) {
			if (numbers[j] < numbers[i]) {
				i = j;
			}
		}

		return i;
	}

	private static boolean isProgramCorrect() {
		boolean isSequentialMinimumSearchCorrect = isMethodCorrect(sequentialMinimumSearch);
		boolean isRecursiveMinimumSearchCorrect = isMethodCorrect(recursiveMinimumSearch);
		boolean isParallelMinimumSearchCorrect = isMethodCorrect(parallelMinimumSearch);

		return isSequentialMinimumSearchCorrect && isRecursiveMinimumSearchCorrect && isParallelMinimumSearchCorrect;
	}

	private static boolean isMethodCorrect(Function<int[], Integer> minimumSearchFunction) {
		int[] randomNumbers = random.ints(100_000_000, 0, 999_999).toArray();

		final int expectedIndex = randomNumbers.length - 1;
		randomNumbers[expectedIndex] = -1;
		int actualIndex = minimumSearchFunction.apply(randomNumbers);

		return actualIndex == expectedIndex;
	}

}
