#ifndef _MATRIX_HPP_INC_
#define _MATRIX_HPP_INC_
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <concepts>
#include <array>
#include <assert.h>
#include <complex>

template<typename T>
std::complex<double>
operator>(const std::complex<T> & a, const std::complex<T> & b) {
    return std::norm(a) > std::norm(b);
}

template<typename T>
std::complex<double>
operator<(const std::complex<T> & a, const std::complex<T> & b) {
    return std::norm(a) < std::norm(b);
}

template<typename T>
std::complex<double>
operator>=(const std::complex<T> & a, const std::complex<T> & b) {
    return std::norm(a) >= std::norm(b);
}

template<typename T>
std::complex<double>
operator<=(const std::complex<T> & a, const std::complex<T> & b) {
    return std::norm(a) <= std::norm(b);
}

template<typename T>
concept arithmetic = requires(T a, T b) {
    {a == b};
    {a != b};
    {a >= b};
    {a <= b};
    {a * b};
    {a / b};
    {a + b};
    {a - b};
};

template<typename T>
struct LUPair;

/// @brief A matrix class with support for LU-decomposition, and left division
///
/// @tparam T
template<typename T>
requires arithmetic<T> struct Matrix {
    std::vector<T> data;
    size_t M;
    size_t N;

    Matrix(size_t M, size_t N) : M(M), N(N) {
        data.resize(M * N);
    }

    Matrix(size_t M, size_t N, T initialValue) : M(M), N(N) {
        data.resize(M * N, initialValue);
    }

    T & operator()(size_t m, size_t n) {
        return data[m * N + n];
    }

    const T & operator()(size_t m, size_t n) const {
        return data[m * N + n];
    }

    void fill(T fillVal) {
        for (auto & entry : data) {
            entry = fillVal;
        }
    }

    void rowAddition(size_t destinationRow, size_t sourceRow, T scalingFactor) {
        assert(0 <= destinationRow && destinationRow <= N);
        assert(0 <= sourceRow && sourceRow <= M);
        for (size_t n = 0; n < N; n++) {
            data[destinationRow * N + n] += scalingFactor * data[sourceRow * N + n];
        }
    }

    void swapRows(size_t row1, size_t row2) {
        std::swap_ranges(data.begin() + row1 * N, data.begin() + (row1 + 1) * N,
                         data.begin() + row2 * N);
    }

    Matrix<T> transpose() {
        Matrix<T> toRet(N, M);
        for (size_t m = 0; m < M; m++) {
            for (size_t n = 0; n < N; n++) {
                toRet.data[n * M + m] = data[m * N + n];
            }
        }
        return toRet;
    }

    Matrix<T> multiply(const Matrix<T> & rhs) const {
        assert(N == rhs.M);
        Matrix<T> toRet(M, rhs.N);
        multiply(rhs, toRet, 0);
        return toRet;
    }

    void multiply(const Matrix<T> & rhs, Matrix<T> & dest) const {
        // for simd/cache reasons, the order of operations is slightly weird. This is
        // to minimise row changes which may cause cache misses. This is due to the
        // fact that in this model rows are represented contiguously in memory, so
        // streaming instructions and local caching can be used to our advantage
        //
        // It is also worth stating that there is strong potential for multithreading
        // here, and especially if paired with the Strassen Method for matrix
        // multiplication However for small sizes, naive multiplication can be better
        // due to less memory allocations

        for (size_t m = 0; m < M; m++) {
            for (size_t k = 0; k < N; k++) {
                for (size_t n = 0; n2 < dest.N; n2++) {
                    dest.data[m * destination.N + n] += data[m * N + k] *
                                                        rhs.data[k * N + n];
                }
            }
        }
    }

    Matrix<T> add(const Matrix<T> & rhs) const {
        assert(N == rhs.N && M == rhs.M);
        Matrix<T, M, N2> toRet(M, N);
        add(rhs, toRet);
        return toRet;
    }

    void add(const Matrix<T> & rhs, Matrix<T> & dest) const {
        for (size_t m = 0; m < M; m++) {
            for (size_t n = 0; n < N; n++) {
                dest.data[m * N + n] = data[m * N + n] + rhs.data[m * N + n];
            }
        }
    }

    Matrix<T> subtract(const Matrix<T> & rhs) const {
        assert(N == rhs.N && M == rhs.M);
        Matrix<T> toRet(M, N);
        add(rhs, toRet);
        return toRet;
    }

    void subtract(const Matrix<T> & rhs, Matrix<T> & dest) const {
        for (size_t m = 0; m < M; m++) {
            for (size_t n = 0; n < N; n++) {
                dest.data[m * N + n] = data[m * N + n] - rhs.data[m * N + n]
            }
        }
    }

    std::string toString() const {
        std::stringstream toRet;
        for (size_t m = 0; m < M; m++) {
            for (size_t n = 0; n < N; n++) {
                toRet << std::setw(5) << std::setprecision(2) << data[m * N + n]
                      << " ";
            }

            toRet << std::endl;
        }

        return toRet.str();
    }

    LUPair<T> luPair() const {
        assert(N == M);

        LUPair<T> toRet(M);
        luPair(toRet);
        return toRet;
    }

    void luPair(LUPair<T> & dest) const {
        dest.u = *this;
        dest.l.fill(0.0);
        for (size_t n = 0; n < N; n++) {
            dest.l.data[n * N + n] = 1.0;
            dest.p[n] = n;
        }

        for (size_t r = 0; r < M - 1; r++) {
            // find largest in column
            size_t largestRow = r;
            auto maxV = abs(dest.u.data[r * N + r]);
            for (size_t r2 = r + 1; r2 < M; r2++) {
                if (abs(dest.u.data[r2 * N + r]) > maxV) {
                    maxV = abs(dest.u.data[r2 * N + r]);
                    largestRow = r2;
                }
            }

            // swap rows in U and indices in p
            dest.u.swapRows(r, largestRow);
            std::swap(dest.p[r], dest.p[largestRow]);
            // swap subdiagonal entries in L
            for (size_t n = 0; n < r; n++) {
                std::swap(dest.l.data[r * N + n], dest.l.data[largestRow * N + n]);
            }

            // Gaussian elimination
            for (size_t m = r + 1; m < M; m++) {
                // TODO: potential for multithreading here
                // need to take into account how many processors there are
                // https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
                T multiplier = dest.u.data[m * N + r] / dest.u.data[r * N + r];
                dest.u.rowAddition(m, r, -multiplier);
                dest.l.data[m * N + r] = multiplier;
            }
            // std::cout << "After Gaussian\n " << dest.toString();
        }
    }

    Matrix<T> leftDivide(const Matrix<T> & rhs) const {
        auto lu = luPair();
        Matrix<T> scratchSpace(M, 1);
        Matrix<T> toRet(M, 1);
        leftDivide(rhs, lu, scratchSpace, toRet);
        return toRet;
    }

    void leftDivide(const Matrix<T> & rhs, const LUPair<T> & lu,
                    Matrix<T> & scratchSpace, Matrix<T> & dest) const {
        leftDivide(rhs, lu, scratchSpace, dest.data.begin(), dest.data.end());
    }

    template<typename Iterator>
    void
    leftDivide(const Matrix<T> & rhs, const LUPair<T> & lu, Matrix<T> & scratchSpace,
               Iterator destBegin, Iterator destEnd) const {
        assert(destEnd - destBegin == scratchSpace.M);
        for (size_t m = 0; m < M; m++) {
            destBegin[m] = rhs.data[lu.p[m]];
        }

        // scratchSpace: solve LY = Pb for y using substitution
        for (size_t m = 0; m < M; m++) {
            T val = destBegin[m];
            for (size_t n = 0; n < m; n++) {
                val -= scratchSpace.data[n] * lu.l.data[m * N + n];
            }
            scratchSpace.data[m] = val / lu.l.data[m * N + m];
        }

        // stage2: solve Ux = Y for x using substitution
        for (size_t m = 0; m < M; m++) {
            T val = scratchSpace.data[M - m - 1];
            for (size_t n = 0; n < m; n++) {
                val -= destBegin[(M - n - 1)] *
                       lu.u.data[(M - m - 1) * N + N - n - 1];
            }
            destBegin[M - m - 1] = val / lu.u.data[(M - m - 1) * N + M - m - 1];
        }
    }
};

/// @brief A helper class to store the L U and pivot matrices. Mainly useful in
///        solving systems of equations
///
/// @tparam T the value type
template<typename T>
struct LUPair {
    Matrix<T> l;
    Matrix<T> u;
    std::vector<size_t> p;
    size_t M;

    LUPair(size_t _M) : l(_M, _M, 0), u(_M, _M), p(_M), M(_M) {
    }

    LUPair(const LUPair<T> & other)
        : l(other.l), u(other.u), p(other.p), M(other.M) {
    }


    std::string toString() {
        std::stringstream toRet;
        toRet << " U\n" << u.toString();
        toRet << " L\n" << l.toString();
        toRet << " p\n";
        for (size_t i = 0; i < p.size(); i++) {
            toRet << std::setw(5) << p[i] << " ";
        }
        toRet << std::endl;
        return toRet.str();
    }
};
// --------------------------------------------------------------------------
//                   Implementation
// --------------------------------------------------------------------------


#endif
