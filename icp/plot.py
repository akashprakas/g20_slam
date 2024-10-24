import numpy as np
import matplotlib.pyplot as plt

def read_matrix_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    true_data = []
    transformed_data = []
    current_matrix = None

    for line in lines:
        line = line.strip()

        if line == "True Data:":
            current_matrix = true_data
        elif line == "Transformed Data:":
            current_matrix = transformed_data
        elif line and current_matrix is not None:
            # Convert space-separated string to list of floats
            current_matrix.append(list(map(float, line.split())))

    return np.array(true_data), np.array(transformed_data)

def plot_data(true_data, transformed_data):
    # Extract x and y coordinates
    true_x, true_y = true_data[0, :], true_data[1, :]
    transformed_x, transformed_y = transformed_data[0, :], transformed_data[1, :]

    # Plot true data
    plt.scatter(true_x, true_y, color='blue', label='True Data', marker='o')
    plt.plot(true_x, true_y, color='blue', linestyle='--')

    # Plot transformed data
    plt.scatter(transformed_x, transformed_y, color='red', label='Transformed Data', marker='x')
    plt.plot(transformed_x, transformed_y, color='red', linestyle='-')

    # Labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('True Data vs Transformed Data')
    plt.legend()

    # Show plot
    plt.show()

if __name__ == "__main__":
    # Read matrices from file
    true_data, transformed_data = read_matrix_from_file('/home/akash/Documents/code/icp_cmake/matrix_log.txt')

    # Plot the data
    plot_data(true_data, transformed_data)
