#include <iostream>
#include "Maths/DynamicMatrix.hpp"
#include <thread>
#include <complex>
#include <chrono>

[[clang::optnone]] int
main() {
    unsigned int nThreads = std::thread::hardware_concurrency();
    std::cout << nThreads << " concurrent threads are supported.\n";

    std::cout << "2 x 2 Matrix\n" << std::endl;
    Matrix<double> myMat(2, 2);
    // std::cin >> d1 >> d2 >> d3 >> d4;
    myMat(0, 0) = 1.0;
    myMat(0, 1) = 1.0;
    myMat(1, 0) = 1.0;
    myMat(1, 1) = 2.0;
    std::cout << " myMat\n" << myMat.toString();

    auto startTime = std::chrono::high_resolution_clock::now();
    auto lu = myMat.luPair();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto timeTaken = std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
                                                                           startTime)
                         .count();

    std::cout << timeTaken << "ms" << std::endl << lu.toString();

    std::cout << "4 x 4 Matrix\n" << std::endl;
    Matrix<double> myMat4(4, 4);
    Matrix<std::complex<double> > myCMat4(4, 4);
    // std::cin >> d1 >> d2 >> d3 >> d4;
    myMat4(0, 0) = 2.0;
    myMat4(0, 1) = 1.0;
    myMat4(0, 2) = 1.0;
    myMat4(0, 3) = 0.0;
    myMat4(1, 0) = 4.0;
    myMat4(1, 1) = 3.0;
    myMat4(1, 2) = 3.0;
    myMat4(1, 3) = 1.0;
    myMat4(2, 0) = 8.0;
    myMat4(2, 1) = 7.0;
    myMat4(2, 2) = 9.0;
    myMat4(2, 3) = 5.0;
    myMat4(3, 0) = 6.0;
    myMat4(3, 1) = 7.0;
    myMat4(3, 2) = 9.0;
    myMat4(3, 3) = 8.0;

    Matrix<double> myCol(4, 1);
    myCol(0, 0) = 1.0;
    myCol(1, 0) = 0.0;
    myCol(2, 0) = 3.0;
    myCol(3, 0) = 0.0;

    myCMat4(0, 0) = 2.0;
    myCMat4(0, 1) = 1.0;
    myCMat4(0, 2) = 1.0;
    myCMat4(0, 3) = 0.0;
    myCMat4(1, 0) = 4.0;
    myCMat4(1, 1) = 3.0;
    myCMat4(1, 2) = 3.0;
    myCMat4(1, 3) = 1.0;
    myCMat4(2, 0) = 8.0;
    myCMat4(2, 1) = 7.0;
    myCMat4(2, 2) = 9.0;
    myCMat4(2, 3) = 5.0;
    myCMat4(3, 0) = 6.0;
    myCMat4(3, 1) = 7.0;
    myCMat4(3, 2) = 9.0;
    myCMat4(3, 3) = 8.0;

    std::cout << " myMat\n" << myCMat4.toString();

    auto luC4 = myCMat4.luPair();
    auto lu4 = myMat4.luPair();
    std::cout << lu4.toString();
    std::cout << luC4.toString();

    startTime = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < 1000000 /* ( unsigned int )-1*/; i++) {
        lu4 = myMat4.luPair();
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
                                                                      startTime)
                    .count();
    std::cout << 1000000 << " 4 x 4 pairs in " << timeTaken << "ms" << std::endl;

    auto solved = myMat4.leftDivide(myCol);
    std::cout << solved.toString();

    Matrix<double> myMat1000(1000, 1000);
    for (size_t i = 0; i < 1000; i++) {
        myMat1000(i, 1000 - 1 - i) = i + 10;
    }
    startTime = std::chrono::high_resolution_clock::now();
    auto lu1000 = myMat1000.luPair();
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
                                                                      startTime)
                    .count();
    std::cout << "1000 x 1000 LU comp in " << timeTaken << "ms" << std::endl;

    // std::cout << lu1000.toString();

    return 0;
}
