#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>
#include <vector>
#include "matrix.h"


enum importance{
    ALL,
    NEVER,
    IMPORTANT
};

class Logger
{
    void print(Vector& x, int step) const;

protected:
    importance imp = ALL;
    uint step = 1;
    uint counter = 0;

public:
    Logger();
    virtual ~Logger();

    void log(std::string message, bool isImportant = false) const;
    void log(std::string message, double num, bool isImportant = false) const;
    void log(std::string message, Vector& vec, bool isImportant = false) const;

    void logPeriodical(std::string message, bool isImportant = false);
    void logPeriodical(std::string message, double num, bool isImportant = false);
    void logPeriodical(std::string message, Vector& vec, bool isImportant = false);

    void setImportance(importance A);
    void setStep(uint newStep);
};
#endif // LOGGER_H
