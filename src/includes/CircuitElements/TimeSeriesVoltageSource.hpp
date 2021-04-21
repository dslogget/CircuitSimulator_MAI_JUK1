#ifndef _TIMESERIESVOLTAGESOURCE_INC_
#define _TIMESERIESVOLTAGESOURCE_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include <math.h>

/// @brief A time series voltage source model
///
/// @tparam T the value type
template<typename T>
struct TimeSeriesVoltageSource : public Component<T> {
    size_t n1 = 0;
    size_t n2 = 0;
    size_t currentIndex = 0;
    size_t lastTimeSeriesIndex = 0;

    std::vector<T> timeSeries;
    std::vector<T> dataSeries;

    T lerp(size_t lowIndex, T timeVal) const {
        T diffTS = timeSeries[(lowIndex + 1) % timeSeries.size()] -
                   timeSeries[lowIndex];
        T diffTV = timeVal - timeSeries[lowIndex];
        T diffDS = dataSeries[(lowIndex + 1) % timeSeries.size()] -
                   dataSeries[lowIndex];
        return dataSeries[lowIndex] + diffDS * diffTV / diffTS;
    }

    void addDynamicStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep) const {
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;
        size_t currentIndexp = currentIndex - 1;
        size_t timeSeriesIndex = lastTimeSeriesIndex;
        T timeMod = std::fmod(currentSolutionIndex * timestep, timeSeries.back());
        while (timeMod > timeSeries[(timeSeriesIndex + 1) % timeSeries.size()] ||
               (timeSeriesIndex != 0 &&
                timeMod < timeSeries[(timeSeriesIndex - 1) % timeSeries.size()])) {
            timeSeriesIndex = (timeSeriesIndex + 1) % timeSeries.size();
        }

        if (n1) {
            stamp.G(n1p, stamp.sizeG_A + currentIndexp) += 1;
            stamp.G(stamp.sizeG_A + currentIndexp, n1p) += 1;
        }

        if (n2) {
            stamp.G(n2p, stamp.sizeG_A + currentIndexp) += -1;
            stamp.G(stamp.sizeG_A + currentIndexp, n2p) += -1;
        }

        stamp.s(stamp.sizeG_A + currentIndexp, 0) += lerp(timeSeriesIndex, timeMod);
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
        while (std::fmod(currentSolutionIndex * timestep, timeSeries.back()) >
                   timeSeries[(lastTimeSeriesIndex + 1) % timeSeries.size()] ||
               (lastTimeSeriesIndex != 0 &&
                std::fmod(currentSolutionIndex * timestep, timeSeries.back()) <
                    timeSeries[(lastTimeSeriesIndex - 1) % timeSeries.size()])) {
            lastTimeSeriesIndex = (lastTimeSeriesIndex + 1) % timeSeries.size();
        }
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        addDynamicStampTo(stamp, solutionVector, 0, 0);
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        // n1 n2 timescale file
        std::regex timeSeriesVoltageSourceRegex = generateRegex("VT", "n n w s",
                                                                true, false);

        TimeSeriesVoltageSource<T> voltageSource;
        std::smatch matches;

        std::regex_search(line, matches, timeSeriesVoltageSourceRegex);

        voltageSource.designator = "VT";
        voltageSource.designator += matches.str(1);

        voltageSource.n1 = std::stoi(matches.str(2));
        voltageSource.n2 = std::stoi(matches.str(3));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));

        T timescale = 0;
        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
            timescale = std::stod(matches.str(4));
        } else {
            static_assert("Unsupported Type");
        }


        voltageSource.currentIndex = ++numCurrents;

        std::regex endOfLineRegex(R"(^(.*)$)");
        std::regex_search(matches.suffix().first, line.cend(), matches,
                          endOfLineRegex);

        std::string seriesPath = matches.str(1);
        voltageSource.readInTimeSeries(timescale, seriesPath);

        elements.dynamicElements.emplace_back(
            std::make_shared<TimeSeriesVoltageSource<T> >(voltageSource));
        elements.nodeComponentMap.insert(
            {{voltageSource.n1, elements.dynamicElements.back()},
             {voltageSource.n2, elements.dynamicElements.back()}});
    }

    void readInTimeSeries(T timescale, const std::string & seriesPath) {
        std::ifstream file(seriesPath);

        std::string line;
        T time;
        T val;

        while (!file.eof()) {
            if (!std::isdigit(file.peek())) {
                std::getline(file, line);
                continue;
            }

            file >> time;

            while (!std::isdigit(file.peek()) &&
                   !(file.peek() == '-' || file.peek() == '+')) {
                file.ignore(1);
            }

            file >> val;

            timeSeries.emplace_back(time * timescale);
            dataSeries.emplace_back(val);
        }
    }
};

#endif
