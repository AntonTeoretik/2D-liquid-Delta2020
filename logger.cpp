#include "logger.h"

using namespace std;

void Logger::setImportance(importance A)
{
    imp = A;
}

void Logger::setStep(uint newStep)
{
    step = newStep;
}

Logger::Logger() {}

Logger::~Logger() {}

void Logger::log(string message, bool isImportant) const
{
    if (imp == NEVER or (imp == IMPORTANT and not isImportant) )
        return;
    std::cout << message << std::endl;
}

void Logger::log(string message, double num, bool isImportant) const
{
    if (imp == NEVER or (imp == IMPORTANT and not isImportant) )
        return;
    std::cout << message << " " << num << std::endl;
}

void Logger::log(string message, Vector &vec, bool isImportant) const
{
    if (imp == NEVER or (imp == IMPORTANT and not isImportant) )
        return;


    std::cout << message << " " << std::endl;
    for(size_t i = 0; i < vec.size(); i++)
        std::cout << vec.getValue(i) << " ";
    std::cout << std::endl;
}

void Logger::logPeriodical(string message, bool isImportant)
{
    counter = (counter >= step) ? 1 : counter + 1;

    if (counter == 1)
        log(message, isImportant);
}

void Logger::logPeriodical(string message, double num, bool isImportant)
{
    counter = (counter >= step) ? 1 : counter + 1;

    if (counter == 1)
        log(message, num, isImportant);
}

void Logger::logPeriodical(string message, Vector &vec, bool isImportant)
{
    counter = (counter >= step) ? 1 : counter + 1;

    if (counter == 1)
        log(message, vec, isImportant);
}
