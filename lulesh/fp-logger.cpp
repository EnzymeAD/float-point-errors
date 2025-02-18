#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "fp-logger.hpp"

constexpr double EPS = 0.0;

class ValueInfo {
public:
  double minRes = std::numeric_limits<double>::max();
  double maxRes = std::numeric_limits<double>::lowest();
  std::vector<double> minOperands;
  std::vector<double> maxOperands;
  unsigned executions = 0;

  double runningSumLog = 0.0;
  unsigned runningCountNonZero = 0;
  double runningSumArith = 0.0;
  unsigned validCount = 0;

  void update(double res, const double *operands, unsigned numOperands) {
    ++executions;

    if (minOperands.empty()) {
      minOperands.resize(numOperands, std::numeric_limits<double>::max());
      maxOperands.resize(numOperands, std::numeric_limits<double>::lowest());
    }
    for (unsigned i = 0; i < numOperands; ++i) {
      if (!std::isnan(operands[i])) {
        minOperands[i] = std::min(minOperands[i], operands[i]);
        maxOperands[i] = std::max(maxOperands[i], operands[i]);
      }
    }

    if (!std::isnan(res)) {
      minRes = std::min(minRes, res);
      maxRes = std::max(maxRes, res);

      double absRes = std::fabs(res);
      runningSumArith += absRes;
      ++validCount;

      if (EPS != 0.0) {
        runningSumLog += std::log(absRes + EPS);
      } else {
        if (absRes != 0.0) {
          runningSumLog += std::log(absRes);
          ++runningCountNonZero;
        }
      }
    }
  }

  double getGeoMean() const {
    if (validCount == 0)
      return 0.0;

    if (EPS != 0.0) {
      return std::exp(runningSumLog / validCount) - EPS;
    } else {
      if (runningCountNonZero == 0) {
        return 0.0;
      }
      return std::exp(runningSumLog / runningCountNonZero);
    }
  }

  double getArithMean() const {
    if (validCount == 0)
      return 0.0;
    return runningSumArith / validCount;
  }

  double getMaxAbs() const {
    return std::max(std::abs(minRes), std::abs(maxRes));
  }
};

class GradInfo {
public:
  double runningSumLog = 0.0;
  unsigned runningCountNonZero = 0;
  double runningSumArith = 0.0;
  unsigned validCount = 0;

  double minGrad = std::numeric_limits<double>::max();
  double maxGrad = std::numeric_limits<double>::lowest();

  void update(double grad) {
    if (!std::isnan(grad)) {
      minGrad = std::min(minGrad, grad);
      maxGrad = std::max(maxGrad, grad);

      double absGrad = std::fabs(grad);

      runningSumArith += absGrad;
      ++validCount;

      if (EPS != 0.0) {
        runningSumLog += std::log(absGrad + EPS);
      } else {
        if (absGrad != 0.0) {
          runningSumLog += std::log(absGrad);
          ++runningCountNonZero;
        }
      }
    }
  }

  double getGeoMean() const {
    if (validCount == 0)
      return 0.0;

    if (EPS != 0.0) {
      return std::exp(runningSumLog / validCount) - EPS;
    } else {
      if (runningCountNonZero == 0) {
        return 0.0;
      }
      return std::exp(runningSumLog / runningCountNonZero);
    }
  }

  double getArithMean() const {
    if (validCount == 0)
      return 0.0;
    return runningSumArith / validCount;
  }

  double getMaxAbs() const {
    return std::max(std::abs(minGrad), std::abs(maxGrad));
  }
};

class ErrorInfo {
public:
  double minErr = std::numeric_limits<double>::max();
  double maxErr = std::numeric_limits<double>::lowest();

  void update(double err) {
    if (!std::isnan(err)) {
      minErr = std::min(minErr, err);
      maxErr = std::max(maxErr, err);
    }
  }
};

class Logger {
private:
  std::unordered_map<std::string, ValueInfo> valueInfo;
  std::unordered_map<std::string, ErrorInfo> errorInfo;
  std::unordered_map<std::string, GradInfo> gradInfo;

public:
  void updateValue(const std::string &id, double res, unsigned numOperands,
                   const double *operands) {
    auto &info = valueInfo.emplace(id, ValueInfo()).first->second;
    info.update(res, operands, numOperands);
  }

  void updateError(const std::string &id, double err) {
    auto &info = errorInfo.emplace(id, ErrorInfo()).first->second;
    info.update(err);
  }

  void updateGrad(const std::string &id, double grad) {
    auto &info = gradInfo.emplace(id, GradInfo()).first->second;
    info.update(grad);
  }

  void print() const {
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (const auto &pair : valueInfo) {
      const auto &id = pair.first;
      const auto &info = pair.second;
      std::cout << "Value:" << id << "\n";
      std::cout << "\tMinRes = " << info.minRes << "\n";
      std::cout << "\tMaxRes = " << info.maxRes << "\n";
      std::cout << "\tExecutions = " << info.executions << "\n";
      std::cout << "\tGeoMeanAbs = " << info.getGeoMean() << "\n";
      std::cout << "\tArithMeanAbs = " << info.getArithMean() << "\n";
      std::cout << "\tMaxAbs = " << info.getMaxAbs() << "\n";
      for (unsigned i = 0; i < info.minOperands.size(); ++i) {
        std::cout << "\tOperand[" << i << "] = [" << info.minOperands[i] << ", "
                  << info.maxOperands[i] << "]\n";
      }
    }

    for (const auto &pair : errorInfo) {
      const auto &id = pair.first;
      const auto &info = pair.second;
      std::cout << "Error:" << id << "\n";
      std::cout << "\tMinErr = " << info.minErr << "\n";
      std::cout << "\tMaxErr = " << info.maxErr << "\n";
    }

    for (const auto &pair : gradInfo) {
      const auto &id = pair.first;
      const auto &info = pair.second;
      std::cout << "Grad:" << id << "\n";
      std::cout << "\tGeoMeanAbs = " << info.getGeoMean() << "\n";
      std::cout << "\tArithMeanAbs = " << info.getArithMean() << "\n";
      std::cout << "\tMaxAbs = " << info.getMaxAbs() << "\n";
    }
  }
};

Logger *logger = nullptr;

void initializeLogger() { logger = new Logger(); }

void destroyLogger() {
  delete logger;
  logger = nullptr;
}

void printLogger() { logger->print(); }

void enzymeLogError(const char *id, double err) {
  assert(logger && "Logger is not initialized");
  logger->updateError(id, err);
}

void enzymeLogGrad(const char *id, double grad) {
  assert(logger && "Logger is not initialized");
  logger->updateGrad(id, grad);
}

void enzymeLogValue(const char *id, double res, unsigned numOperands,
                    double *operands) {
  assert(logger && "Logger is not initialized");
  logger->updateValue(id, res, numOperands, operands);
}
