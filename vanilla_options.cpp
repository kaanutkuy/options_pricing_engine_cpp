#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include "option_class.hpp"
#include "utility_functions.hpp"

using namespace std;


class EuropeanOption : public Option {
public:
    // Constructor
    EuropeanOption(double currentVal, double strike, double vol, double divRate, double rfRate, double time, bool isCall)
        : Option(currentVal, strike, vol, divRate, rfRate, time, isCall) {}

    double blackScholes_europeanOption() {
        double d1 = UtilityFunctions::d1_(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility, dividendRate);
        double d2 = UtilityFunctions::d2_(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility, dividendRate);

        if (isCallOption) {
            optionPrice = exp(-riskFreeRate * timeToExpiry) * (underlyingPrice * exp((riskFreeRate - dividendRate) * timeToExpiry) * UtilityFunctions::cdf(d1) - strikePrice * UtilityFunctions::cdf(d2));
            return optionPrice;
        }
        else if (!isCallOption) {
            optionPrice = exp(-riskFreeRate * timeToExpiry) * (strikePrice * UtilityFunctions::cdf(-d2) - underlyingPrice * exp((riskFreeRate - dividendRate) * timeToExpiry) * UtilityFunctions::cdf(-d1));
            return optionPrice;
        }
        return 0;
    }

    double monteCarlo_europeanOption(int timeSteps, int simulations) {
        double dt = timeToExpiry / timeSteps;
        double nu = (riskFreeRate - dividendRate - 0.5 * volatility * volatility); 
        double erddt = exp((riskFreeRate - dividendRate) * dt);
        double egamma = exp((2 * (riskFreeRate - dividendRate) + volatility * volatility) * dt) - 2 * erddt + 1;
        double beta1 = -1;              // beta coefficient on delta control variate
        double beta2 = -0.5;            // beta coefficient on gamma control variate
        double option_sum = 0.0;

        // uses mersenne twister (pseudo-random) and one seed
        // NOTE: might try to have more homogenous randomness using quasi-random generators (like sobol)
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);  // normal distribution with mean 0 & variance 1

        for (int j = 0; j < simulations; j++) {
            double St1 = underlyingPrice;
            double St2 = underlyingPrice;
            double cv1 = 0.0;
            double cv2 = 0.0;

            for (int i = 0; i < timeSteps; i++) {
                double number = distribution(gen);  // get random number from normal distribution

                double delta1 = UtilityFunctions::delta(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double delta2 = UtilityFunctions::delta(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma1 = UtilityFunctions::gamma(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma2 = UtilityFunctions::gamma(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);

                double Stn1 = St1 * exp(nu*dt + volatility*sqrt(dt)*number);
                double Stn2 = St2 * exp(nu*dt + volatility*sqrt(dt)*(-number));

                cv1 = cv1 + delta1 * (Stn1 - St1 * erddt) + 
                            delta2 * (Stn2 - St2 * erddt);
                cv2 = cv2 + gamma1 * ((Stn1 - St1)*(Stn1 - St1) - St1*St1*egamma) +
                            gamma2 * ((Stn2 - St2)*(Stn2 - St2) - St2*St2*egamma);

                St1 = Stn1;
                St2 = Stn2;
            }
            double option_T = 0.5 * (max(0.0, ((isCallOption ? 1 : -1) * (St1 - strikePrice))) + beta1*cv1 + \
                                    max(0.0, ((isCallOption ? 1 : -1) * (St2 - strikePrice))) + beta2*cv2);  // take the average of the option prices 
            option_sum += option_T;
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    }

    double monteCarlo_europeanOption_heston(int timeSteps, int simulations, double kappa, double theta, double sigma, double rho) {
        double dt = timeToExpiry / timeSteps;
        double erddt = exp((riskFreeRate - dividendRate) * dt);
        double beta1 = -1;          // beta coefficient on delta control variate
        double beta2 = -0.5;        // beta coefficient on gamma control variate
        double option_sum = 0.0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0, 1.0);  // normal distribution with mean 0 & variance 1

        for (int j = 0; j < simulations; j++) {
            double St1 = underlyingPrice;
            double St2 = underlyingPrice;
            double vt1 = volatility;
            double vt2 = volatility;
            double cv1 = 0.0;
            double cv2 = 0.0;

            for (int i = 0; i < timeSteps; i++) {
                double zS = distribution(gen);  // get random number from normal distribution
                double zV = distribution(gen);

                double delta1 = UtilityFunctions::delta(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double delta2 = UtilityFunctions::delta(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma1 = UtilityFunctions::gamma(St1, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);
                double gamma2 = UtilityFunctions::gamma(St2, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry - i*dt, isCallOption);

                double tempS1 = St1; // store previous St values before they get updated in hestonDynamics
                double tempS2 = St2;

                UtilityFunctions::hestonDynamics(St1, vt1, riskFreeRate, dividendRate, kappa, theta, sigma, rho, dt, zS, zV);
                UtilityFunctions::hestonDynamics(St2, vt2, riskFreeRate, dividendRate, kappa, theta, sigma, rho, dt, -zS, -zV); // antithetic variate

                cv1 += delta1 * (St1 - tempS1 * erddt) + 
                       delta2 * (St2 - tempS2 * erddt);
                cv2 += gamma1 * ((St1 - tempS1)*(St1 - tempS1) - St1*St1*(exp((2*(riskFreeRate - dividendRate) + vt1*vt1) * dt) - 2*erddt + 1)) +
                       gamma2 * ((St2 - tempS2)*(St2 - tempS2) - St2*St2*(exp((2*(riskFreeRate - dividendRate) + vt2*vt2) * dt) - 2*erddt + 1));
            }
            double option_T = 0.5 * (max(0.0, ((isCallOption ? 1 : -1) * (St1 - strikePrice))) + beta1*cv1 + \
                                     max(0.0, ((isCallOption ? 1 : -1) * (St2 - strikePrice))) + beta2*cv2);
            option_sum += option_T;
        }
        double option_value = exp(-riskFreeRate * timeToExpiry) * (option_sum / simulations);
        optionPrice = option_value;
        return option_value;
    };

    // function to calculate and store the Greeks. If it is not called, the Greeks will be set to 0 for every new instance of the class
    void calculateGreeks() {
        delta = UtilityFunctions::delta(underlyingPrice, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry, isCallOption);
        gamma = UtilityFunctions::gamma(underlyingPrice, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry, isCallOption);
        vega = UtilityFunctions::vega(underlyingPrice, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry, isCallOption);
        rho = UtilityFunctions::rho(underlyingPrice, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry, isCallOption);
        theta = UtilityFunctions::theta(underlyingPrice, strikePrice, riskFreeRate, dividendRate, volatility, timeToExpiry, isCallOption);
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Option Delta is: ";
        cout << getDelta() << endl;
        cout << "Option Gamma is: ";
        cout << getGamma() << endl;
        cout << "Option Vega is: ";
        cout << getVega() << endl;
        cout << "Option Rho is: ";
        cout << getRho() << endl;
        cout << "Option Theta is: ";
        cout << getTheta() << endl;
    }
};

class AmericanOption : public Option {
protected:
    double timeStepsNumber = 100;   // placeholder until set in one of the pricer methods

public:
    // Constructor
    AmericanOption(double currentVal, double strike, double vol, double divRate, double rfRate, double time, bool isCall)
        : Option(currentVal, strike, vol, divRate, rfRate, time, isCall) {}

    double crr_americanOption(int timeSteps) {
        timeStepsNumber = timeSteps;
        double dt = timeToExpiry / timeSteps;

        double nu = (riskFreeRate - dividendRate - 0.5 * volatility * volatility);
        double lambda = sqrt(3/2);  // Kamrad and Ritchken best convergence rate sqrt(3/2)
        double pu = 1 / (2 * lambda * lambda) + 0.5 * (nu / (lambda * volatility)) * sqrt(dt);
        double pd = 1 / (2 * lambda * lambda) - 0.5 * (nu / (lambda * volatility)) * sqrt(dt);
        double pm = 1 - pu - pd;

        double u = exp(lambda * volatility * sqrt(dt));
        double p = 1/u;

        double S[timeSteps+1][2*timeSteps+1];
        double C[timeSteps+1][2*timeSteps+1];
        // vector<vector<double> > S(N+1, vector<double>(2*N+1));
        // vector<vector<double> > C(N+1, vector<double>(2*N+1));

        for (int i = timeSteps; i >= 0; i--) {
            for (int j = -i; j <= i; j++) {
                S[i][j] = underlyingPrice * pow(u, j);
            }
        }

        for (int j = timeSteps; j >= -timeSteps; j--) {
            C[timeSteps][j] = max(0.0, ((isCallOption ? 1 : -1) * (S[timeSteps][j] - strikePrice)));
        }

        for (int i = timeSteps-1; i >= 0; i--) {
            for (int j = i; j >= -i; j--) {
                C[i][j] = max((exp(-riskFreeRate * dt) * (pu*C[i+1][j+1] + pm*C[i+1][j] + pd*C[i+1][j-1])), ((isCallOption ? 1 : -1) * (S[i][j] - strikePrice)));
            }
        }
        optionPrice = C[0][0];
        return C[0][0];
    };

    double ls_americanOption(int timeSteps) {
        if (timeSteps % 2 == 0) {  // use odd time-steps for correct convergence 
            timeSteps += 1;
        }
        timeStepsNumber = timeSteps;
        double dt = timeToExpiry / timeSteps;

        double d1 = UtilityFunctions::d1_(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility, dividendRate);
        double d2 = UtilityFunctions::d2_(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility, dividendRate);

        double pu = UtilityFunctions::peizer_pratt_cdf(d2, timeSteps);
        double pd = UtilityFunctions::peizer_pratt_cdf(d1, timeSteps);

        double u = exp((riskFreeRate - dividendRate) * dt) * (pd / pu);
        double d = exp((riskFreeRate - dividendRate) * dt) * ((1 - pd) / (1 - pu));

        vector<vector<double> > S(timeSteps+1, vector<double>(timeSteps+1));
        vector<vector<double> > C(timeSteps+1, vector<double>(timeSteps+1));

        for (int i = 0; i <= timeSteps; i++) {
            for (int j = 0; j <= i; j++) {
                S[i][j] = underlyingPrice * pow(u, j) * pow(d, i-j);
            }
        }

        for (int j = 0; j <= timeSteps; j++) {
            C[timeSteps][j] = max(0.0, ((isCallOption? 1 : -1) * (S[timeSteps][j] - strikePrice)));
        }

        for (int i = timeSteps-1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                C[i][j] = max((exp(-riskFreeRate * dt) * (pu*C[i+1][j+1] + (1-pu)*C[i+1][j])), ((isCallOption ? 1 : -1) * (S[i][j] - strikePrice)));
            }
        }
        optionPrice = C[0][0];

        delta = (C[1][1] - C[1][0]) / (S[1][1] - S[1][0]);
        gamma = (((C[2][2] - C[2][1]) / (S[2][2] - S[2][1])) - ((C[2][1] - C[2][0]) / (S[2][1] - S[2][0]))) / (0.5 * (S[2][2] - S[2][0]));

        return C[0][0];
    };

    // function to calculate and store the Greeks. If it is not called, the Greeks will be set to 0 for every new instance of the class
    void calculateGreeks() {
        double dV = 0.0001 * volatility;
        double dR = 0.0001 * riskFreeRate;
        double dT = 0.0001 * timeToExpiry;

        double v_C1 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility + dV, dividendRate, timeStepsNumber, isCallOption);
        double v_C2 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate, volatility - dV, dividendRate, timeStepsNumber, isCallOption);
        vega = (v_C1 - v_C2) / (2 * dV);

        double r_C1 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate + dR, volatility, dividendRate, timeStepsNumber, isCallOption);
        double r_C2 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry, riskFreeRate - dR, volatility, dividendRate, timeStepsNumber, isCallOption);
        rho = (r_C1 - r_C2) / (2 * dR);

        double t_C1 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry + dT, riskFreeRate, volatility, dividendRate, timeStepsNumber, isCallOption);
        double t_C2 = UtilityFunctions::lr_americanOption_forGreeks(underlyingPrice, strikePrice, timeToExpiry - dT, riskFreeRate, volatility, dividendRate, timeStepsNumber, isCallOption);
        theta = (t_C1 - t_C2) / (2 * dT);
    }

    void displayInfo() override {
        cout << "Option Price is: ";
        cout << getOptionPrice() << endl;
        cout << "Option Delta is: ";
        cout << getDelta() << endl;
        cout << "Option Gamma is: ";
        cout << getGamma() << endl;
        cout << "Option Vega is: ";
        cout << getVega() << endl;
        cout << "Option Rho is: ";
        cout << getRho() << endl;
        cout << "Option Theta is: ";
        cout << getTheta() << endl;
    }
};
