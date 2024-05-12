#include <iostream>
#include <vector>
#include <cmath>
#include <random>
// #include <stdexcept>

#include "utility_functions.hpp"

using namespace std;


// use this for european IF volatility is LOW (implied volatility is low)
double blackScholes_europeanOption(double S, double K, double T, double r, double sig, double dividend, bool callOption) {
    double d1 = UtilityFunctions::d1_(S, K, T, r, sig, dividend);
    double d2 = UtilityFunctions::d2_(S, K, T, r, sig, dividend);

    if (callOption) {
        return exp(-dividend * T) * S * UtilityFunctions::cdf(d1) - K * exp(-r * T) * UtilityFunctions::cdf(d2);
    }
    else if (!callOption) {
        return exp(-dividend * T) * S * (UtilityFunctions::cdf(d1) - 1) - K * exp(-r * T) * (UtilityFunctions::cdf(d2) - 1);
    }
    else { // LOOK INTO ERROR CATCHING
        //throw std::invalid_argument("Invalid option type. Please specify either CALL (true) or PUT (false).");
        cout << "Invalid option type. Please input 'true' for CALL or 'false' PUT." << endl;
        return 0;
    }
}


int main() {
    double sig = 0.2;  // volatility 
    double S = 100;  // current stock price
    double dividend = 0.03;  // dividend
    double K = 100;  // option strike price
    double r = 0.06;  // risk-free interest rate
    double T = 1;  // time until expiration (in years)
    bool callOption = true;  // true for Call, false for Put

    cout << "European " << (callOption ? "Call" : "Put") << " Option Value: " << blackScholes_europeanOption(S, K, T, r, sig, dividend, callOption) << endl;

    return EXIT_SUCCESS;
}