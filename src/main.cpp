
#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

using Matrix = algebra::Matrix;

int main(int argc, char **argv)
{
    if (0) // make false to run unit-tests
    {
        // debug section
        Matrix matrix{{-1, 1.5, -1.75, -2}, {-2, 2.5, -2.75, -3}, {3, 3.5, -3.75, -4}, {4, 4.5, 4.75, -5}};
        Matrix inverse{algebra::inverse(matrix)};
        //EXPECT_NEAR(inverse[0][0], 0.16, 0.03);
        //EXPECT_NEAR(inverse[1][1], 3.31, 0.03);
        //EXPECT_NEAR(inverse[3][1], 2.67, 0.03);
        //EXPECT_NEAR(inverse[0][3], 0, 0.03);
        for (int i = 0; i < matrix.size(); i ++) 
        {
            for (int j = 0; j < matrix.size(); j ++)
            {
                std::cout << inverse[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
    
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}