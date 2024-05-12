#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include "option_class.hpp"
#include "utility_functions.hpp"

using namespace std;


class AsianOptionFixed : public Option {
public:
    // Constructor
    AsianOptionFixed(double currentVal, double strike, double vol, double divRate, double rfRate, double time, bool isCall)
        : Option(currentVal, strike, vol, divRate, rfRate, time, isCall) {}

    double monteCarlo_asianOption_geometric(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double product = 1.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);  // get random number from normal distribution
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                product *= St;
            }
            double geometric_avg = pow(product, 1.0/timeSteps);

            option_sum += max(0.0, ((isCallOption ? 1 : -1) * (geometric_avg - strikePrice)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    };

    double monteCarlo_asianOption_arithmetic(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double sum = 0.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);  // get random number from normal distribution
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                sum += St;
            }
            double arithmetic_avg = sum / timeSteps;

            option_sum += max(0.0, ((isCallOption ? 1 : -1) * (arithmetic_avg - strikePrice)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    };

    double monteCarlo_asianOption_geometricControlVariate(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double sum = 0.0;
            double product = 1.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                sum += St;
                product *= St;
            }
            double arithmetic_avg = sum / timeSteps;
            double geometric_average = pow(product, 1.0/timeSteps);

            option_sum += max(0.0, ((isCallOption ? 1 : -1)*(arithmetic_avg - strikePrice) - (isCallOption ? 1 : -1)*(geometric_average - strikePrice)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations) + this->monteCarlo_asianOption_geometric(timeSteps, simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};

// Asian Option with floating strike
class AsianOptionFloating : public FloatingOption {
public:
    // Constructor
    AsianOptionFloating(double currentVal, double vol, double divRate, double rfRate, double time, bool isCall)
        : FloatingOption(currentVal, vol, divRate, rfRate, time, isCall) {}

    double monteCarlo_asianOption_geometric(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double product = 1.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);  // get random number from normal distribution
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                product *= St;
            }
            double geometric_avg = pow(product, 1.0/timeSteps);

            option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - geometric_avg)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    };

    double monteCarlo_asianOption_arithmetic(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double sum = 0.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);  // get random number from normal distribution
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                sum += St;
            }
            double arithmetic_avg = sum / timeSteps;

            option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - arithmetic_avg)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    };

    double monteCarlo_asianOption_geometricControlVariate(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double sum = 0.0;
            double product = 1.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                sum += St;
                product *= St;
            }
            double arithmetic_avg = sum / timeSteps;
            double geometric_average = pow(product, 1.0/timeSteps);

            option_sum += max(0.0, ((isCallOption ? 1 : -1)*(St - arithmetic_avg) - (isCallOption ? 1 : -1)*(St - geometric_average)));
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations) + this->monteCarlo_asianOption_geometric(timeSteps, simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};


class LookbackOptionFixed : public Option {
public:
    // Constructor
    LookbackOptionFixed(double currentVal, double strike, double vol, double divRate, double rfRate, double time, bool isCall)
        : Option(currentVal, strike, vol, divRate, rfRate, time, isCall) {}

    double monteCarlo_lookbackOption_fixed(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5 * volatility * volatility);
        double erddt = exp((riskFreeRate - dividendRate) * dt);
        double egamma = exp((2 * (riskFreeRate - dividendRate) + volatility * volatility) * dt) - 2 * erddt + 1;
        double beta1 = -1;
        double beta2 = -0.5;
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St1 = underlyingPrice;
            double St2 = underlyingPrice;
            double cv1 = 0.0;
            double cv2 = 0.0;
            double St1_min = St1;
            double St1_max = St1;
            double St2_min = St2;
            double St2_max = St2;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);

                double delta1 = UtilityFunctions::delta(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double delta2 = UtilityFunctions::delta(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma1 = UtilityFunctions::gamma(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma2 = UtilityFunctions::gamma(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);

                double Stn1 = St1 * exp(nu*dt + volatility * sqrt(dt) * number);
                double Stn2 = St2 * exp(nu*dt + volatility * sqrt(dt) * (-number));

                cv1 = cv1 + delta1 * (Stn1 - St1 * erddt) + 
                            delta2 * (Stn2 - St2 * erddt);
                cv2 = cv2 + gamma1 * ((Stn1 - St1)*(Stn1 - St1) - St1*St1*egamma) +
                            gamma2 * ((Stn2 - St2)*(Stn2 - St2) - St2*St2*egamma);

                St1 = Stn1;
                St2 = Stn2;
                St1_min = min(St1, St1_min);
                St1_max = max(St1, St1_max);
                St2_min = min(St2, St2_min);
                St2_max = max(St2, St2_max);
            }
            if (isCallOption) {
                double option_T = 0.5 * (max(0.0, St1_max - strikePrice) + beta1*cv1 + \
                                        max(0.0, St2_max - strikePrice) + beta2*cv2);  // take the average of the option prices 
                option_sum += option_T;
            }
            else if (!isCallOption) {
                double option_T = 0.5 * (max(0.0, strikePrice - St1_min) + beta1*cv1 + \
                                        max(0.0, strikePrice - St2_min) + beta2*cv2);  // take the average of the option prices 
                option_sum += option_T;
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};

class LookbackOptionFloating : public FloatingOption {
public:
    // Constructor
    LookbackOptionFloating(double currentVal, double vol, double divRate, double rfRate, double time, bool isCall)
        : FloatingOption(currentVal, vol, divRate, rfRate, time, isCall) {}

    double monteCarlo_lookbackOption_floating(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            double St_min = St;
            double St_max = St;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);
                St_min = min(St, St_min);
                St_max = max(St, St_max);
            }
            if (isCallOption) {
                option_sum += St - St_min;
            }
            else if (!isCallOption) {
                option_sum += St_max - St;
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};


class BarrierOption : public Option {
protected:
    double barrierPrice;    // Price of the set barrier

public:
    // Constructor
    BarrierOption(double currentVal, double strike, double vol, double divRate, double rfRate, double time, double barrier, bool isCall)
        : Option(currentVal, strike, vol, divRate, rfRate, time, isCall), barrierPrice(barrier) {}

    // Getter methods
    double getBarrierPrice() const { return barrierPrice; }

    // Setter methods 
    void setBarrierPrice(double barrier) { barrierPrice = barrier; }
};

class BarrierDownIn : public BarrierOption {
public:
    // Constructor
    BarrierDownIn(double currentVal, double strike, double vol, double divRate, double rfRate, double time, double barrier, bool isCall)
        : BarrierOption(currentVal, strike, vol, divRate, rfRate, time, barrier, isCall) {}

    double monteCarlo_barrierOption_downIn(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            bool crossed = false;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);

                if (St <= barrierPrice) {
                    crossed = true;
                }
            }
            if (crossed) {
                option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - strikePrice)));
            }
            else {
                option_sum += 0;
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};

class BarrierUpIn : public BarrierOption {
public:
    // Constructor
    BarrierUpIn(double currentVal, double strike, double vol, double divRate, double rfRate, double time, double barrier, bool isCall)
        : BarrierOption(currentVal, strike, vol, divRate, rfRate, time, barrier, isCall) {}

    double monteCarlo_barrierOption_upIn(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            bool crossed = false;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);

                if (St >= barrierPrice) {
                    crossed = true;
                }
            }
            if (crossed) {
                option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - strikePrice)));
            }
            else {
                option_sum += 0;
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};

class BarrierDownOut : public BarrierOption {
public:
    // Constructor
    BarrierDownOut(double currentVal, double strike, double vol, double divRate, double rfRate, double time, double barrier, bool isCall)
        : BarrierOption(currentVal, strike, vol, divRate, rfRate, time, barrier, isCall) {}

    double monteCarlo_barrierOption_downOut(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            bool crossed = false;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);

                if (St <= barrierPrice) {
                    crossed = true;
                }
            }
            if (crossed) {
                option_sum += 0;
            }
            else {
                option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - strikePrice)));
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};

class BarrierUpOut : public BarrierOption {
public:
    // Constructor
    BarrierUpOut(double currentVal, double strike, double vol, double divRate, double rfRate, double time, double barrier, bool isCall)
        : BarrierOption(currentVal, strike, vol, divRate, rfRate, time, barrier, isCall) {}

    double monteCarlo_barrierOption_upOut(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5*volatility*volatility);
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (int j = 0; j < simulations; j++) {
            double St = underlyingPrice;
            bool crossed = false;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);
                St *= exp(nu*dt + volatility*sqrt(dt)*number);

                if (St >= barrierPrice) {
                    crossed = true;
                }
            }
            if (crossed) {
                option_sum += 0;
            }
            else {
                option_sum += max(0.0, ((isCallOption ? 1 : -1) * (St - strikePrice)));
            }
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Greeks information currently not available for exotic options, sorry!" << endl;
    }
};
