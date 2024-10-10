import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.util.*;

public class Graph {
	/**
	 * Write a method that takes in as input a graph [matrix] and
	 * returns a number of different values (integer) depending
	 * on whether the matrix is valid / invalid, as
	 * defined/described above. Note that the matrix will be
	 * represented as a 2- dimensional array of integers
	 *
	 * @param arr The 2-dimensional array that needs to be tested.
	 * @return 0 if invalid, 1 if valid.
	 */
	public static int isValidMatrix(int[][] arr) {

		//Check the matrix is square.
		int rows = arr.length;
		if (rows == 0) {
			return 0;
		}

		int columns = arr[0].length;
		if (columns == 0) {
			return 0;
		}

		if (rows != columns) {
			return 0;
		}

		//Check the matrix is symmetrical
		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < columns; x++) {
				int a = arr[y][x]; //rows columns
				int b = arr[x][y]; //column row

				//Check A is either 1 or 0
				if (a != 1 && a != 0) {
					return 0;
				}

				//Are the values symmetrical?
				if (a != b) {
					return 0;
				}

				//Check that the values are 0, diagonally
				if (y == x) {
					if (a != 0) {
						return 0;
					}
				}
			}
		}

		return 1;
	}

	/**
	 * Write a method that creates the initial starting colouring
	 * arrangement (the colours MUST be numbered 1, 2, 3 or
	 * 4). Note that the method should take in N as a parameter
	 * and return the arrangement as output. The starting
	 * arrangement should be in vector format e.g. ArrayList
	 * of integers
	 *
	 * @param N Size of row & column
	 * @return Vector List of nodes with their colours set linearly.
	 */
	public static List<Integer> initialStartingPoint(int N) {
		if (N <= 0) throw new RuntimeException("N must be bigger than 0");

		int size = N;

		//Create new Arraylist
		ArrayList<Integer> arrayList = new ArrayList<>(size);

		//Apply algorithm to assign values.
		for (int i = 0; i < size; i++) {
			arrayList.add(rand(1, 4));
		}

		return arrayList;
	}

	/**
	 * Write a method that changes the value of a single random
	 * element of the arrangement to another random value,
	 * ensuring that the new value is a valid colour. The method
	 * should take the colouring arrangement as input and
	 * returns the changed colouring arrangement as output (a
	 * copy MUST be returned). Note that the arrangement will
	 * be a vector format, such as an ArrayList of integers
	 * @param node
	 * @param matrix
	 * @param solution
	 * @return
	 */
	public static int colourClashesOfaNode(int node, int[][] matrix, List<Integer> solution) {
		if (isValidMatrix(matrix) != 1) {
			throw new RuntimeException("Matrix is not valid");
		}

		if (isSolutionValid(solution, matrix.length) != 1) {
			throw new RuntimeException("Solution is invalid");
		}

		int currentColour = solution.get(node);
		int clash = 0;

		for (int k = 0; k < matrix.length; k++) {
			int current = matrix[node][k];

			//If it is 0, it does not matter.
			if (current != 1 || node == k) {
				continue;
			}

			int thisColour = solution.get(k);

			if (currentColour == thisColour) {
				clash++;
			}
		}

		return clash;
	}

	/**
	 * Write a method that takes in a graph [matrix] and a solution
	 * to the Graph Colouring Problem and returns the fitness
	 * value based on the description of the fitness function
	 * above
	 */
	public static int fitnessFunction(int[][] matrix, List<Integer> solution) {
		int fitness = 0;
		if (isValidMatrix(matrix) != 1) {
			throw new RuntimeException("Matrix is invalid");
		}

		if (isSolutionValid(solution, matrix.length) != 1) {
			throw new RuntimeException("Solution is invalid");
		}

		for (int i = 0; i < matrix.length; i++) {
			fitness += colourClashesOfaNode(i, matrix, solution);
		}

		return fitness;
	}

	/**
	 * Write a method that changes the value of a single random
	 * element of the arrangement to another random value,
	 * ensuring that the new value is a valid colour. The method
	 * should take the colouring arrangement as input and
	 * returns the changed colouring arrangement as output (a
	 * copy MUST be returned). Note that the arrangement will
	 * be a vector format, such as an ArrayList of integers.
	 *
	 * @return
	 */
	public static List<Integer> smallChangeOperator(List<Integer> solution) {
		int index = rand(0, solution.size() - 1);
		int colour = rand(1, 4);
		int currentColour = solution.get(index);

		if (colour == currentColour){
			return smallChangeOperator(solution);
		}

		List<Integer> newList = new ArrayList<>(solution);
		newList.set(index, colour);

		return newList;
	}

	/** //GA Solution Finder
	 * Write a class that when given a graph [matrix] and the
	 * number of iterations it applies a single-population
	 * heuristic search algorithm of your choice and returns
	 * the colouring arrangement as output. This class can
	 * contain the components from the previous questions
	 * @param iterations
	 * @param matrix
	 * @return
	 */
	public static List<Integer> solvingTheProblem(int iterations, int[][] matrix) {
		int singlePopulation = iterations / 10; //What is the size of the single population.
		int generations = iterations - singlePopulation; //How many generations?
		boolean enforceIttr = true; //Force the iteration limit? if false, the limit can be exceeded.
		boolean alwaysRun = false; //Ignore the generation limit?

		if (isValidMatrix(matrix) != 1) {
			throw new RuntimeException("Matrix is invalid");
		}

		if (iterations < 1) {
			throw new RuntimeException("the fitness function calls (iterations) being greater than one");
		}

		Population[] pop = new Population[singlePopulation];

		for (int i = 0; i < singlePopulation; i++) {
			Population current = new Population(matrix);
			pop[i] = current;
		}

		for (int a = 0; a < generations || alwaysRun; a++) {
			try {
				sort(pop, enforceIttr, iterations);
			} catch (IllegalAccessException e) {
				System.out.println("Ran out of iterations!");
				return pop[0].getSolution();
			}

			if (pop[0].getFitness() == 0) {
				return pop[0].getSolution();
			}

			for (int i = 0; i < pop.length - 1; i++) {
				int breed = rand(1, 2);

				if (breed != 2) {
					continue;
				}

				Population current = pop[i];

				Population next = pop[rand(0, matrix.length/10)+1];

				int midpoint = current.getSolution().size() / 2;

				int[] currentGenetics = current.getSolution().stream().mapToInt(Integer::intValue).toArray();
				int[] nextGenetics = next.getSolution().stream().mapToInt(Integer::intValue).toArray();

				for (int x = midpoint; x < currentGenetics.length; x++) {
					int temp = currentGenetics[x];
					currentGenetics[x] = nextGenetics[x];
					nextGenetics[x] = temp;
				}

				current.setSolution(Arrays.stream(currentGenetics).boxed().toList());
			}

			pop[pop.length - 1].setSolution(new ArrayList<>(pop[0].getSolution()));

			for (int i = 1; i < pop.length; i++) {
				Population current = pop[i];
				current.setSolution(smallChangeOperator(current.getSolution()));
				current.setSolution(smallChangeOperator(current.getSolution()));
			}
		}

		return pop[0].getSolution();
	}

	//RMHC Solution Finder
	public static List<Integer> solvingTheProblem2(int iterations, int[][] matrix) {
		Population population = new Population(matrix);
		population.setSolution(initialStartingPoint(matrix.length));

		int bestFitness = -1;
		List<Integer> bestSolution = null;
		for (int i = 0; i < iterations; i++) {
			population.setSolution(smallChangeOperator(population.getSolution()));
			population.setSolution(smallChangeOperator(population.getSolution()));
			population.setSolution(smallChangeOperator(population.getSolution()));

			population.updateFitness();

			if (bestSolution == null){
				bestSolution = new ArrayList<>(population.getSolution());
				bestFitness = population.getFitness();
				continue;
			}

			if (population.getFitness() > bestFitness){
				bestSolution = new ArrayList<>(population.getSolution());
				bestFitness = population.getFitness();

			}
		}

		return population.getSolution();
	}

	public static void sort(Population[] populations, boolean enforceMax, int max) throws IllegalAccessException {
		for (Population population : populations) {

			if (enforceMax){
				if (Population.calls+1 >= max){
					throw new IllegalAccessException("Too many calls");
				}
			}

			population.updateFitness();
		}

		for (int x = 0; x < populations.length; x++) {
			boolean swapped = false;
			for (int i = 1; i < populations.length; i++) {
				Population current = populations[i - 1];
				Population next = populations[i];

				//if current is smaller than next
				if (current.getFitness() > next.getFitness()) {
					populations[i - 1] = next;
					populations[i] = current;
					swapped = true;
				}
			}

			if (!swapped) {
				break;
			}
		}
	}

	public static int rand(int aa, int bb) {
		return Util.UI(aa, bb);
	}

	public static int isSolutionValid(List<Integer> solution, int N) {
		int size = N;

		if (solution.size() != size) {
			return 0;
		}

		return 1;
	}

	private static class Population {
		private List<Integer> solution;
		private int fitness = 0;

		private int[][] matrix;

		private static int calls;

		public Population(int[][] matrix) {
			this.matrix = matrix;
			this.setSolution(initialStartingPoint(this.matrix.length));
		}

		public List<Integer> getSolution() {
			return solution;
		}

		public void setSolution(List<Integer> solution) {
			this.solution = solution;
		}

		public int getFitness() {
			return fitness;
		}

		public void updateFitness(){
			fitness = fitnessFunction(this.matrix, this.getSolution());
			calls++;
		}

		public static int getCalls() {
			return calls;
		}
	}

	static class test {
		public static void main(String[] args) {
			System.out.println("---------------[testMatrix]---------------");
			testMatrix();

			System.out.println("---------------[testStartingPoint]---------------");
			testStartingPoint();

			System.out.println("---------------[testFitnessFunction]---------------");
			testFitnessFunction();
			System.out.println("---------------------------------------------------");

			int[][] matrix2 = {
					{0,	1,	1,	0,	1,	0,	0,	1,	0,	1},
					{1,	0,	0,	1,	0,	0,	0,	0,	0,	0},
					{1,	0,	0,	0,	0,	0,	1,	0,	0,	0},
					{0,	1,	0,	0,	1,	0,	0,	0,	0,	1},
					{1,	0,	0,	1,	0,	0,	0,	1,	0,	1},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	1,	0,	0,	0,	0,	1,	0,	0},
					{1,	0,	0,	0,	1,	0,	1,	0,	0,	1},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	1},
					{1,	0,	0,	1,	1,	0,	0,	1,	1,	0}};

			List<Integer> solution = solvingTheProblem2(100, matrix2);

			System.out.println("Fitness: " + fitnessFunction(matrix2, solution) + " Top: " + solution);
			System.out.println("calls: " + Population.calls);
		}

		public static void testFitnessFunction() {
			int[][] matrix = {{0, 1, 0, 0, 1}, {1, 0, 1, 0, 1}, {0, 1, 0, 1, 0}, {0, 0, 1, 0, 1}, {1, 1, 0, 1, 0,}};
			int[][] matrix1 = {{0, 0, 0, 0, 1}, {0, 0, 1, 0, 0}, {0, 1, 0, 1, 0}, {0, 0, 1, 0, 0}, {1, 0, 0, 0, 0,}};
			ArrayList<Integer> solution = initialStartingPoint(5);

			System.out.print("Test Case 1: \t\t\t");
			System.out.println(fitnessFunction(matrix, solution) == 2 ? "Passed" : "Failed");

			System.out.print("Test Case 2: \t\t\t");
			System.out.println(fitnessFunction(matrix1, solution) == 2 ? "Passed" : "Failed");

			int[][] matrix2 = {
					{0, 1, 0, 0, 1},
					{1, 0, 1, 0, 1},
					{0, 1, 0, 1, 0},
					{0, 0, 1, 0, 1},
					{1, 1, 0, 1, 0,}};

			int[] sol2 = new int[]{1, 3, 1, 3, 1};

			System.out.print("Test Case 3: \t\t\t");
			System.out.println(fitnessFunction(matrix2, Arrays.stream(sol2).boxed().toList()) == 2 ? "Passed" : "Failed");
		}

		public static void testMatrix() {

			//Non-Square
			int[][] invalidMatrix1 = {{0, 1, 1}, {1, 0, 1}};

			System.out.print("Non-Square: \t\t\t");
			System.out.println(Graph.isValidMatrix(invalidMatrix1) == 0 ? "Passed" : "Fail");

			//Non-Symmetrical
			int[][] invalidMatrix2 = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};

			System.out.print("Non-Symmetrical: \t\t");
			System.out.println(Graph.isValidMatrix(invalidMatrix2) == 0 ? "Passed" : "Fail");

			//Non-Binary Values
			int[][] invalidMatrix3 = {{0, 2, 1}, {2, 0, 3}, {1, 3, 0}};

			System.out.print("Non-Binary Values: \t\t");
			System.out.println(Graph.isValidMatrix(invalidMatrix3) == 0 ? "Passed" : "Fail");

			//Diagonal Non-Zeros
			int[][] invalidMatrix4 = {{1, 0, 1}, {0, 1, 0}, {1, 0, 1}};

			System.out.print("Diagonal Non-Zeros: \t");
			System.out.println(Graph.isValidMatrix(invalidMatrix4) == 0 ? "Passed" : "Fail");

			int[][] validMatrix1 = {{0, 1, 0, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}, {1, 0, 1, 0}};

			int[][] validMatrix2 = {{0, 1, 1, 0}, {1, 0, 1, 1}, {1, 1, 0, 1}, {0, 1, 1, 0}};

			int[][] validMatrix3 = {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}};

			System.out.print("Valid Matrix 1: \t\t");
			System.out.println(Graph.isValidMatrix(validMatrix1) == 1 ? "Passed" : "Fail");

			System.out.print("Valid Matrix 2: \t\t");
			System.out.println(Graph.isValidMatrix(validMatrix2) == 1 ? "Passed" : "Fail");

			System.out.print("Valid Matrix 3: \t\t");
			System.out.println(Graph.isValidMatrix(validMatrix3) == 1 ? "Passed" : "Fail");
		}

		public static void testStartingPoint() {

			try {
				// Test case 3: N = 0 (Should return an error arrangement)

				List<Integer> testCase3 = Graph.initialStartingPoint(0);

				System.out.print("N0 Test Case 3: \t\t");
				System.out.println("Fail");
			} catch (RuntimeException e) {
				System.out.print("N0 Test Case 3: \t\t");
				System.out.println("Passed");
			}
		}

		public static ArrayList<Integer> initialStartingPoint(int N) {
			if (N <= 0) throw new RuntimeException("N must be bigger than 0");

			int size = N;

			//Create new Arraylist
			ArrayList<Integer> arrayList = new ArrayList<>(size);

			//Apply algorithm to assign values.
			for (int i = 0; i < size; i++) {
				arrayList.add((i % 4) + 1);
			}

			return arrayList;
		}
	}


	public class Util {
		//Shared random object
		static private Random rand;

		//Create a uniformly distributed random integer between aa and bb inclusive
		static public int UI(int aa, int bb) {
			int a = Math.min(aa, bb);
			int b = Math.max(aa, bb);
			if (rand == null) {
				rand = new Random();
				rand.setSeed(System.nanoTime());
			}
			int d = b - a + 1;
			int x = rand.nextInt(d) + a;
			return (x);
		}

		//Create a uniformly distributed random double between a and b inclusive
		static public double UR(double a, double b) {
			if (rand == null) {
				rand = new Random();
				rand.setSeed(System.nanoTime());
			}
			return ((b - a) * rand.nextDouble() + a);
		}

		static public ArrayList<Double> ReadNumberFile(String filename) {
			ArrayList<Double> res = new ArrayList<Double>();
			Reader r;
			try {
				r = new BufferedReader(new FileReader(filename));
				StreamTokenizer stok = new StreamTokenizer(r);
				stok.parseNumbers();
				stok.nextToken();
				while (stok.ttype != StreamTokenizer.TT_EOF) {
					if (stok.ttype == StreamTokenizer.TT_NUMBER) {
						res.add(stok.nval);
					}
					stok.nextToken();
				}
			} catch (Exception E) {
				System.out.println("+++ReadFile: " + E.getMessage());
			}
			return (res);
		}
	}
}
