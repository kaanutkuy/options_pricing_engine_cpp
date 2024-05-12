#ifndef OPTION_CLASS_H
#define OPTION_CLASS_H

#include <iostream>
#include <cmath>

#include "utility_functions.hpp"

using namespace std;

class Option {
protected:
    double underlyingPrice;     // Current value of the underlying asset
    double strikePrice;         // Strike price of the option
    double volatility;          // Volatility of the underlying asset
    double dividendRate;        // Dividend rate of the underlying asset
    double riskFreeRate;        // Risk-free interest rate
    double timeToExpiry;        // Time until expiry of the option
    bool isCallOption;          // Flag indicating if it's a call option (true) or put option (false)

    double delta=0.0;           // Sensitivity of option price to change in the price of underlying (dC/dS)
    double gamma=0.0;           // Sensitivity of option price to a change in delta (change in the change in underlying price) (d^2C/dS^2)
    double vega=0.0;            // Sensitivity of option price to change in volatility (dC/dv)
    double rho=0.0;             // Sensitivity of option price to change in time (dC/dt)
    double theta=0.0;           // Sensitivity of option price to change in risk-free interest rate (dC/dr)

    double optionPrice=0.0;     // Value of the option. Set to 0 until priced by appropriate methods

public:
    // Constructor
    Option(double currentVal, double strike, double vol, double div, double r, double time, bool isCall)
        : underlyingPrice(currentVal), strikePrice(strike), volatility(vol), dividendRate(div),
          riskFreeRate(r), timeToExpiry(time), isCallOption(isCall) {}

    // Getter methods
    double getUnderlyingPrice() const { return underlyingPrice; }
    double getStrikePrice() const { return strikePrice; }
    double getVolatility() const { return volatility; }
    double getDividendRate() const { return dividendRate; }
    double getRiskFreeRate() const { return riskFreeRate; }
    double getTimeToExpiry() const { return timeToExpiry; }
    bool isCall() const { return isCallOption; }
    double getDelta() const { return delta; }
    double getGamma() const { return gamma; }
    double getVega() const { return vega; }
    double getRho() const { return rho; }
    double getTheta() const { return theta; }
    double getOptionPrice() const { return optionPrice; }

    // Setter methods (if needed)
    void setUnderlyingPrice(double val) { underlyingPrice = val; }
    void setStrikePrice(double strike) { strikePrice = strike; }
    void setVolatility(double vol) { volatility = vol; }
    void setDividendRate(double divRate) { dividendRate = divRate; }
    void setRiskFreeRate(double rfRate) { riskFreeRate = rfRate; }
    void setTimeToExpiry(double time) { timeToExpiry = time; }
    void setOptionType(bool isCall) { isCallOption = isCall; }

    virtual void displayInfo() {}
};

// Option class for exotic options where the strike price is not predetermined
class FloatingOption {
protected:
    double underlyingPrice;     // Current value of the underlying asset
    double volatility;          // Volatility of the underlying asset
    double dividendRate;        // Dividend rate of the underlying asset
    double riskFreeRate;        // Risk-free interest rate
    double timeToExpiry;        // Time until expiry of the option
    bool isCallOption;          // Flag indicating if it's a call option (true) or put option (false)

    double delta=0.0;           // Sensitivity of option price to change in the price of underlying (dC/dS)
    double gamma=0.0;           // Sensitivity of option price to a change in delta (change in the change in underlying price) (d^2C/dS^2)
    double vega=0.0;            // Sensitivity of option price to change in volatility (dC/dv)
    double rho=0.0;             // Sensitivity of option price to change in time (dC/dt)
    double theta=0.0;           // Sensitivity of option price to change in risk-free interest rate (dC/dr)

    double optionPrice=0.0;     // Value of the option. Set to 0 until priced by appropriate methods

public:
    // Constructor
    FloatingOption(double currentVal, double vol, double div, double r, double time, bool isCall)
        : underlyingPrice(currentVal), volatility(vol), dividendRate(div),
          riskFreeRate(r), timeToExpiry(time), isCallOption(isCall) {}

    // Getter methods
    double getUnderlyingPrice() const { return underlyingPrice; }
    double getVolatility() const { return volatility; }
    double getDividendRate() const { return dividendRate; }
    double getRiskFreeRate() const { return riskFreeRate; }
    double getTimeToExpiry() const { return timeToExpiry; }
    bool isCall() const { return isCallOption; }
    double getDelta() const { return delta; }
    double getGamma() const { return gamma; }
    double getVega() const { return vega; }
    double getRho() const { return rho; }
    double getTheta() const { return theta; }
    double getOptionPrice() const { return optionPrice; }

    // Setter methods (if needed)
    void setUnderlyingPrice(double val) { underlyingPrice = val; }
    void setVolatility(double vol) { volatility = vol; }
    void setDividendRate(double divRate) { dividendRate = divRate; }
    void setRiskFreeRate(double rfRate) { riskFreeRate = rfRate; }
    void setTimeToExpiry(double time) { timeToExpiry = time; }
    void setOptionType(bool isCall) { isCallOption = isCall; }

    virtual void displayInfo() {}
};
#endif