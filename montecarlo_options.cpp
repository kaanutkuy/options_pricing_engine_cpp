#include <iostream>
#include <vector>
#include <cmath>
#include <random>
// #include <boost/random/sobol.hpp>
// #include <boost/random/random_device.hpp>
#include "utility_functions.hpp"

using namespace std;


// use if implied volatility is HIGH (else, use black scholes)
// monte carlo with antithetic & delta-gamma based variance reduction
double mc_euoropeanOption_naive(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double erddt = exp((r - dividend)*dt);
    double egamma = exp((2*(r - dividend) + sig*sig) * dt) - 2*erddt + 1;
    double beta1 = -1; 
    double beta2 = -0.5;
    double option_sum = 0.0;

    // uses mersenne twister (pseudo-random) and one seed, might try to have better
    // randomness using quasi-random generators (like sobol)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);  // normal distribution with mean 0 & variance 1

    for (int j = 0; j < simulations; j++) {  // for each simulation
        double St1 = S;
        double St2 = S;
        double cv1 = 0.0;
        double cv2 = 0.0;

        for (int i = 0; i < N; i++) {  // for each time step
            double number = distribution(gen);  // get random number from normal distribution

            double delta1 = UtilityFunctions::delta(St1, K, r, dividend, sig, T - i*dt, callOption);
            double delta2 = UtilityFunctions::delta(St2, K, r, dividend, sig, T - i*dt, callOption);
            double gamma1 = UtilityFunctions::gamma(St1, K, r, dividend, sig, T - i*dt, callOption);
            double gamma2 = UtilityFunctions::gamma(St2, K, r, dividend, sig, T - i*dt, callOption);

            double Stn1 = St1 * exp(nu*dt + sig*sqrt(dt)*number);
            double Stn2 = St2 * exp(nu*dt + sig*sqrt(dt)*(-number));

            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt) + 
                        delta2 * (Stn2 - St2 * erddt);
            cv2 = cv2 + gamma1 * ((Stn1 - St1)*(Stn1 - St1) - St1*St1*egamma) +
                        gamma2 * ((Stn2 - St2)*(Stn2 - St2) - St2*St2*egamma);

            St1 = Stn1;
            St2 = Stn2;
        }
        double option_T = 0.5 * (max(0.0, ((callOption ? 1 : -1) * (St1 - K))) + beta1*cv1 + \
                                 max(0.0, ((callOption ? 1 : -1) * (St2 - K))) + beta2*cv2);  // take the average of the option prices 
        option_sum += option_T;
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

// ASIAN OPTIONS
double mc_asianOption_geometric(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        double product = 1.0;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);  // get random number from normal distribution
            St *= exp(nu*dt + sig*sqrt(dt)*number);
            product *= St;
        }
        double geometric_avg = pow(product, 1.0/N);

        option_sum += max(0.0, ((callOption ? 1 : -1) * (geometric_avg - K)));
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_asianOption_arithmetic(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        double sum = 0.0;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);  // get random number from normal distribution
            St *= exp(nu*dt + sig*sqrt(dt)*number);
            sum += St;
        }
        double arithmetic_avg = sum / N;

        option_sum += max(0.0, ((callOption ? 1 : -1) * (arithmetic_avg - K)));
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_asianOption_geometricControlVariate(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        double sum = 0.0;
        double product = 1.0;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);
            sum += St;
            product *= St;
        }
        double arithmetic_avg = sum / N;
        double geometric_average = pow(product, 1.0/N);

        option_sum += max(0.0, ((callOption ? 1 : -1)*(arithmetic_avg - K) - (callOption ? 1 : -1)*(geometric_average - K)));
    }
    double option_value = exp(-r * T) * (option_sum / simulations) + mc_asianOption_geometric(S, K, T, r, sig, dividend, N, simulations, callOption);
    return option_value;
}

// LOOKBACK OPTIONS
double mc_lookbackOption_fixedStrike(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double erddt = exp((r - dividend)*dt);
    double egamma = exp((2*(r - dividend) + sig*sig) * dt) - 2*erddt + 1;
    double beta1 = -1; 
    double beta2 = -0.5;
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St1 = S;
        double St2 = S;
        double cv1 = 0.0;
        double cv2 = 0.0;
        double St1_min = St1;
        double St1_max = St1;
        double St2_min = St2;
        double St2_max = St2;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            double delta1 = UtilityFunctions::delta(St1, K, r, dividend, sig, T - i*dt, callOption);
            double delta2 = UtilityFunctions::delta(St2, K, r, dividend, sig, T - i*dt, callOption);
            double gamma1 = UtilityFunctions::gamma(St1, K, r, dividend, sig, T - i*dt, callOption);
            double gamma2 = UtilityFunctions::gamma(St2, K, r, dividend, sig, T - i*dt, callOption);

            double Stn1 = St1 * exp(nu*dt + sig*sqrt(dt)*number);
            double Stn2 = St2 * exp(nu*dt + sig*sqrt(dt)*(-number));

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
        if (callOption) {
            double option_T = 0.5 * (max(0.0, St1_max - K) + beta1*cv1 + \
                                 max(0.0, St2_max - K) + beta2*cv2);  // take the average of the option prices 
            option_sum += option_T;
        }
        else if (!callOption) {
            double option_T = 0.5 * (max(0.0, K - St1_min) + beta1*cv1 + \
                                 max(0.0, K - St2_min) + beta2*cv2);  // take the average of the option prices 
            option_sum += option_T;
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_lookbackOption_floatingStrike(double S, double T, double r, double sig, double dividend, int N, int simulations, bool callOption) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        double St_min = St;
        double St_max = St;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);
            St_min = min(St, St_min);
            St_max = max(St, St_max);
        }
        if (callOption) {
            option_sum += St - St_min;
        }
        else if (!callOption) {
            option_sum += St_max - St;
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

// BARRIER OPTIONS
double mc_barrierOption_downOut(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption, double barrier) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        bool crossed = false;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);

            if (St <= barrier) {
                crossed = true;
            }
        }
        if (crossed) {
            option_sum += 0;
        }
        else {
            option_sum += max(0.0, ((callOption ? 1 : -1) * (St - K)));
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_barrierOption_upOut(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption, double barrier) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        bool crossed = false;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);

            if (St >= barrier) {
                crossed = true;
            }
        }
        if (crossed) {
            option_sum += 0;
        }
        else {
            option_sum += max(0.0, ((callOption ? 1 : -1) * (St - K)));
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_barrierOption_downIn(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption, double barrier) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        bool crossed = false;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);

            if (St <= barrier) {
                crossed = true;
            }
        }
        if (crossed) {
            option_sum += max(0.0, ((callOption ? 1 : -1) * (St - K)));
        }
        else {
            option_sum += 0;
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

double mc_barrierOption_upIn(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption, double barrier) {
    double dt = T/N;
    double nu = (r - dividend - 0.5*sig*sig);
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int j = 0; j < simulations; j++) {
        double St = S;
        bool crossed = false;

        for (int i = 0; i < N; i++) {
            double number = distribution(gen);
            St *= exp(nu*dt + sig*sqrt(dt)*number);

            if (St >= barrier) {
                crossed = true;
            }
        }
        if (crossed) {
            option_sum += max(0.0, ((callOption ? 1 : -1) * (St - K)));
        }
        else {
            option_sum += 0;
        }
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}


// Heston model, updates the volatility and underlying price changes during path iterations
void hestonDynamics(double& St, double& vt, double r, double dividend, double kappa, double theta, double sigma, double rho, double dt, double z1, double z2) {
    double dWt = sqrt(dt) * z1;
    double dZt = sqrt(dt) * (rho * z1 + sqrt(1 - rho * rho) * z2); 
    
    double psi = sigma * sqrt(vt) * dWt;
    
    vt += kappa * (theta - vt) * dt + sigma * sqrt(vt) * dZt;
    St *= exp((r - dividend - 0.5 * vt) * dt + sqrt(vt) * (rho * dWt + sqrt(1 - rho * rho) * dZt));
}

// monte carlo with stochastic volatility HESTON Model
double mc_euoropeanOption_heston(double S, double K, double T, double r, double sig, double dividend, int N, int simulations, bool callOption, double kappa, double theta, double sigma, double rho) {
    double dt = T/N;
    double erddt = exp((r - dividend) * dt);
    double beta1 = -1;      // beta coefficient on delta control variate
    double beta2 = -0.5;    // beta coefficient on gamma control variate
    double option_sum = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);  // normal distribution with mean 0 & variance 1

    for (int j = 0; j < simulations; j++) {  // for each simulation
        double St1 = S;
        double St2 = S;
        double vt1 = sig;
        double vt2 = sig;
        double cv1 = 0.0;
        double cv2 = 0.0;

        for (int i = 0; i < N; i++) {  // for each time step
            double zS = distribution(gen);  // get random number from normal distribution
            double zV = distribution(gen);

            double delta1 = UtilityFunctions::delta(St1, K, r, dividend, vt1, T - i*dt, callOption);
            double delta2 = UtilityFunctions::delta(St2, K, r, dividend, vt2, T - i*dt, callOption);
            double gamma1 = UtilityFunctions::gamma(St1, K, r, dividend, vt1, T - i*dt, callOption);
            double gamma2 = UtilityFunctions::gamma(St2, K, r, dividend, vt2, T - i*dt, callOption);

            double tempS1 = St1; // store previous St values before they get updated in hestonDynamics
            double tempS2 = St2;

            hestonDynamics(St1, vt1, r, dividend, kappa, theta, sigma, rho, dt, zS, zV);
            hestonDynamics(St2, vt2, r, dividend, kappa, theta, sigma, rho, dt, -zS, -zV); // antithetic variate

            cv1 += delta1 * (St1 - tempS1 * erddt) + 
                   delta2 * (St2 - tempS2 * erddt);
            cv2 += gamma1 * ((St1 - tempS1)*(St1 - tempS1) - St1*St1*(exp((2*(r - dividend) + vt1*vt1) * dt) - 2*erddt + 1)) +
                   gamma2 * ((St2 - tempS2)*(St2 - tempS2) - St2*St2*(exp((2*(r - dividend) + vt2*vt2) * dt) - 2*erddt + 1));
        }
        double option_T = 0.5 * (max(0.0, ((callOption ? 1 : -1) * (St1 - K))) + beta1*cv1 + \
                                 max(0.0, ((callOption ? 1 : -1) * (St2 - K))) + beta2*cv2);  // take the average of the option prices 
        option_sum += option_T;
    }
    double option_value = exp(-r * T) * (option_sum / simulations);
    return option_value;
}

// TODO:
// Monte Carlo options ARCH/GARCH volatilty models


int main() {
    double vol = 0.2;         // volatility
    double dividend = 0.03;   // dividend pay per share
    double S = 100;           // current stock price
    double K = 105;           // option strike price
    double r = 0.06;          // risk-free interest rate
    double T = 1;             // time until expiration (in years)
    int N = 10;               // number of time steps
    int simulations = 10000;  // number of monte carlo simulations
    bool callOption = false;   // true for Call, false for Put

    double barrier = 107;      // set for BARRIER options

    // Heston Model parameters after calibration //
    double kappa = 4.0;       // mean-reversion speed
    double theta = 0.25;      // long-term mean volatility
    double sigma = 0.05;      // volatily of volatility
    double rho = 0.2;         // correlation of stock and volatility Wiener processes

    cout << "European " << (callOption ? "Call" : "Put") << " Option Value: " << mc_euoropeanOption_naive(S, K, T, r, vol, dividend, N, simulations, callOption) << endl;
    cout << "Delta of " << (callOption ? "Call " : "Put ") << UtilityFunctions::delta(S, K, r, dividend, vol, T, callOption) << endl;

    cout << "European " << (callOption ? "Call" : "Put") << " Option Value with Heston Stochastic volatility: " << mc_euoropeanOption_heston(S, K, T, r, vol, dividend, N, simulations, callOption, kappa, theta, sigma, rho) << endl;

    cout << "Geometric Asian " << (callOption ? "Call" : "Put") << " Option Value: " << mc_asianOption_geometric(S, K, T, r, vol, dividend, N, simulations, callOption) << endl;
    cout << "Arithmetic Asian " << (callOption ? "Call" : "Put") << " Option Value: " << mc_asianOption_arithmetic(S, K, T, r, vol, dividend, N, simulations, callOption) << endl;
    cout << "Arithmetic Asian with Geometric Control Variate " << (callOption ? "Call" : "Put") << " Option Value: " << mc_asianOption_geometricControlVariate(S, K, T, r, vol, dividend, N, simulations, callOption) << endl;

    cout << "Lookback Fixed Strike " << (callOption ? "Call" : "Put") << " Option Value: " << mc_lookbackOption_fixedStrike(S, K, T, r, vol, dividend, N, simulations, callOption) << endl;
    cout << "Lookback Floating Strike " << (callOption ? "Call" : "Put") << " Option Value: " << mc_lookbackOption_floatingStrike(S, T, r, vol, dividend, N, simulations, callOption) << endl;

    cout << "Barrier Down and Out " << (callOption ? "Call" : "Put") << " Option Value: " << mc_barrierOption_downOut(S, K, T, r, vol, dividend, N, simulations, callOption, barrier) << endl;
    cout << "Barrier Up and Out " << (callOption ? "Call" : "Put") << " Option Value: " << mc_barrierOption_upOut(S, K, T, r, vol, dividend, N, simulations, callOption, barrier) << endl;
    cout << "Barrier Down and In " << (callOption ? "Call" : "Put") << " Option Value: " << mc_barrierOption_downIn(S, K, T, r, vol, dividend, N, simulations, callOption, barrier) << endl;
    cout << "Barrier Up and In " << (callOption ? "Call" : "Put") << " Option Value: " << mc_barrierOption_upIn(S, K, T, r, vol, dividend, N, simulations, callOption, barrier) << endl;

    return EXIT_SUCCESS;
}
