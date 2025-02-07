#include "hw1.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <random>

namespace algebra {
    using Matrix = std::vector<std::vector<double>>;
    std :: mt19937 gen;
    Matrix zeros(size_t n, size_t m) {
        return Matrix(n, std::vector<double>(m, 0));
    }
    Matrix ones(size_t n, size_t m) {
        return Matrix(n, std::vector<double>(m, 1));
    }
    Matrix random(size_t n, size_t m, double min, double max) {
        if (min > max)
            throw std::logic_error("min cannot be greater than max");
        Matrix matrix(n, std::vector<double>(m));
        for (auto& row : matrix)
            for (auto& elem : row)
                elem = min + gen() * (max - min) / gen.max();
        return matrix;
    }
    void show(const Matrix& matrix) {
        for (const auto& row : matrix) {
            for (const auto& elem : row)
                std::cout << std::setw(3) << elem;
            std::cout << std::endl;
        }
    }
    Matrix multiply(const Matrix& matrix, double c) {
        Matrix result(matrix);
        for (auto& row : result)
            for (auto& elem : row)
                elem *= c;
        return result;
    }
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.size() == 0 || matrix2.size() == 0)
            return Matrix();
        if (matrix1[0].size() != matrix2.size())
            throw std::logic_error("matrices with wrong dimensions cannot be multiplied");
        Matrix result(matrix1.size(), std::vector<double>(matrix2[0].size()));
        for (size_t i = 0; i < result.size(); i++)
            for (size_t j = 0; j < result[i].size(); j++)
                for (size_t k = 0; k < matrix1[0].size(); k++)
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
        return result;
    }
    Matrix sum(const Matrix& matrix, double c) {
        Matrix result(matrix);
        for (auto& row : result)
            for (auto& elem : row)
                elem += c;
        return result;
    }
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty() ^ matrix2.empty())
            throw std::logic_error("sum of an empty matrix");
        if (matrix1.empty() && matrix2.empty())
            return Matrix();
        if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size())
            throw std::logic_error("matrices with wrong dimensions cannot be summed");
        Matrix result(matrix1);
        for (size_t i = 0; i < result.size(); i++)
            for (size_t j = 0; j < result[i].size(); j++)
                result[i][j] += matrix2[i][j];
        return result;
    }
    Matrix transpose(const Matrix& matrix) {
        if (matrix.empty())
            return Matrix();
        Matrix result(matrix[0].size(), std::vector<double>(matrix.size()));
        for (size_t i= 0; i < result.size(); i++)
            for (size_t j = 0; j < result[i].size(); j++)
                result[i][j] = matrix[j][i];
        return result;
    }
    Matrix minor(const Matrix& matrix, size_t row, size_t col) {
        Matrix result(matrix.size() - 1, std::vector<double>(matrix[0].size() - 1));
        for (size_t i = 0; i < matrix.size(); i++)
            for (size_t j = 0; j < matrix[i].size(); j++)
                if (i != row && j != col)
                    result[i < row ? i : i - 1][j < col ? j : j - 1] = matrix[i][j];
        return result;
    }
    double determinant(const Matrix& matrix) {
        if (matrix.empty())
            return 1;
        if (matrix.size() != matrix[0].size())
            throw std::logic_error("non-square matrices have no determinant");
        if (matrix.size() == 1)
            return matrix[0][0];
        double det = 0;
        for (size_t i = 0; i < matrix.size(); i++)
            det += (i % 2 ? -1 : 1) * matrix[0][i] * determinant(minor(matrix, 0, i));
        return det;
    }
    Matrix inverse(const Matrix& matrix) {
        if (matrix.empty())
            return Matrix();
        if (matrix.size() != matrix[0].size())
            throw std::logic_error("non-square matrices have no inverse");
        double det = determinant(matrix);
        if (std::abs(det) < 1e-6)
            throw std::logic_error("non-square matrices have no inverse");
        Matrix result(matrix.size(), std::vector<double>(matrix.size()));
        for (size_t i = 0; i < matrix.size(); i++)
            for (size_t j = 0; j < matrix[i].size(); j++)
                result[i][j] = ((i + j) % 2 ? -1 : 1) * determinant(minor(matrix, j, i)) / det;
        return result;
    }
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis) {
        if (axis == 0) {
            if (matrix1[0].size() != matrix2[0].size())
                throw std::logic_error("matrices with wrong dimensions cannot be concatenated");
            Matrix result(matrix1);
            result.insert(result.end(), matrix2.begin(), matrix2.end());
            return result;
        }
        if (axis == 1) {
            if (matrix1.size() != matrix2.size())
                throw std::logic_error("matrices with wrong dimensions cannot be concatenated");
            Matrix result;
            for (size_t i = 0; i < matrix1.size(); i++) {
                result.push_back(matrix1[i]);
                result.back().insert(result.back().end(), matrix2[i].begin(), matrix2[i].end());
            }
            return result;
        }
        throw std::logic_error("axis must be 0 or 1");
    }
    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
        if (r1 >= matrix.size() || r2 >= matrix.size())
            throw std::logic_error("r1 or r2 inputs are out of range");
        Matrix result(matrix);
        std::swap(result[r1], result[r2]);
        return result;
    }
    Matrix ero_multiply(const Matrix& matrix, size_t r, double c) {
        if (r >= matrix.size())
            throw std::logic_error("r input is out of range");
        Matrix result(matrix);
        for (auto& elem : result[r])
            elem *= c;
        return result;
    }
    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2) {
        if (r1 >= matrix.size() || r2 >= matrix.size() || r1 < 0 || r2 < 0)
            throw std::logic_error("r1 or r2 inputs are out of range");
        Matrix result(matrix);
        for (size_t i = 0; i < result[r1].size(); i++)
            result[r2][i] += c * result[r1][i];
        return result;
    }
    Matrix upper_triangular(const Matrix& matrix) {
        Matrix result(matrix);
        if (matrix.empty())
            return Matrix();
        if (matrix.size() != matrix[0].size())
            throw std::logic_error("non-square matrices have no upper triangular form");
        for (size_t i = 0; i < result.size(); i++) {
            if (result[i][i] == 0) {
                size_t j = i + 1;
                for (; j < result.size(); j++)
                    if (result[j][i] != 0)
                        break;
                if (j == result.size())
                    continue;
                result = ero_swap(result, i, j);
            }
            for (size_t j = i + 1; j < result.size(); j++) {
                double c = -result[j][i] / result[i][i];
                result = ero_sum(result, i, c, j);
            }
        }
        return result;
    }
}