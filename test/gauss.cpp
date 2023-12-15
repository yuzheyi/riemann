#include <iostream>
#include <vector>

// 打印矩阵
void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

// 高斯消元法求解增广矩阵
void gaussianElimination(std::vector<std::vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size() - 1; // 不包括增广列

    for (int i = 0; i < rows; ++i) {
        // 将当前列的第i行元素调整为1
        double divisor = matrix[i][i];
        for (int j = 0; j <= cols; ++j) {
            matrix[i][j] /= divisor;
        }

        // 将其他行的第i列元素消为0
        for (int k = 0; k < rows; ++k) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int j = 0; j <= cols; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }
    }
}

int main() {
    std::vector<std::vector<double>> augmentedMatrix = {
        {1, 0, 0, 1, 0},
        {2, 2, 2, -1, 1},
        {4, 8, 16, 1, 1},
        {8, 24, 72, -1, 2}
    };

    std::cout << "Original Augmented Matrix:" << std::endl;
    printMatrix(augmentedMatrix);

    gaussianElimination(augmentedMatrix);

    std::cout << "\nAugmented Matrix after Gaussian Elimination:" << std::endl;
    printMatrix(augmentedMatrix);

    std::cout << "\nSolution:" << std::endl;
    for (int i = 0; i < augmentedMatrix.size(); ++i) {
        std::cout << "x" << (i + 1) << " = " << augmentedMatrix[i].back() << std::endl;
    }

    return 0;
}
