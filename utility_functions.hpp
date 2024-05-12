#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

#include <cmath> // Include necessary math functions
#include <vector>

using namespace std;


namespace UtilityFunctions {
    
    // probability density function for the standard normal distribution
    double pdf(double x) {
        return (1.0/sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
    }

    // cumulative density function for the standard normal distribution
    double cdf(double x) {
        return 0.5 * (1 + erf(x / sqrt(2.0)));
    }

    double d1_(double S, double K, double T, double r, double sig, double dividend) {
        return (log(S/K) + T*(r - dividend + (0.5 * sig * sig))) / (sig * sqrt(T));
    }

    double d2_(double S, double K, double T, double r, double sig, double dividend) {
        return (log(S/K) + T*(r - dividend - (0.5 * sig * sig))) / (sig * sqrt(T));
    }

    // sensitivity to change in the price of underlying (dC/dS)
    double delta(double S, double K, double r, double q, double sig, double T, bool callOption) {
        double d1 = d1_(S, K, T, r, sig, q);

        if (callOption) {
            return exp(-q * T) * cdf(d1);
        }
        return exp(-q * T) * (cdf(d1) - 1);
    }

    // sensitivity to a change in delta (change in the change in underlying price) (d^2C/dS^2)
    double gamma(double S, double K, double r, double q, double sig, double T, bool callOption) {
        double d1 = d1_(S, K, T, r, sig, q);
        return exp(-q * T) * pdf(d1) / (sig * S * sqrt(T));
    }

    // sensitivity to change in volatility (dC/dv)
    double vega(double S, double K, double r, double q, double sig, double T, bool callOption) {
        double d1 = d1_(S, K, T, r, sig, q);
        return exp(-q * T) * S * sqrt(T) * pdf(d1);
    }

    // sensitivity to change in time (dC/dt)
    double theta(double S, double K, double r, double q, double sig, double T, bool callOption) {
        double d1 = d1_(S, K, T, r, sig, q);
        double d2 = d2_(S, K, T, r, sig, q);

        if (callOption) {
            return (-S * pdf(d1) * sig) / (2 * sqrt(T)) - (r * K * exp(-r * T) * cdf(d2)) + (q * S * exp(-q * T) * cdf(d1));
        }
        return (-S * pdf(d1) * sig) / (2 * sqrt(T)) + (r * K * exp(-r * T) * cdf(-d2)) - (q * S * exp(-q * T) * cdf(-d1));
    }

    // sensitivity to change in risk-free interest rate (dC/dr)
    double rho(double S, double K, double r, double q, double sig, double T, bool callOption) {
        double d2 = d2_(S, K, T, r, sig, q);

        if (callOption) {
            return K * T * exp(-r * T) * cdf(d2);
        }
        return -K * T * exp(-r * T) * cdf(-d2);
    }

    // sign function
    double sign(double number) {
        if (number < 0) {
            return -1;
        }
        else if (number > 0) {
            return 1;
        }
        return 0;
    }

    double peizer_pratt_cdf(double number, double N) {
        double term = pow((number / (N + 1.0/3.0 + 0.1/(N+1))), 2) * (N + 1.0/6.0);
        return 0.5 + 0.5 * sign(number) * sqrt(1 - exp(-term));
    }

    // Heston model, updates the volatility and underlying price changes during path iterations
    // NOTE: Not available currently, needs historical data investigation & the use of discretization
    // and maximum likelihood estimation to estimate the parameters (kappa, theta, sigma, rho)
    void hestonDynamics(double& underlyingCurrent, double& volatilityCurrent, double riskFreeRate, double dividendRate, double kappa, double theta, double sigma, double rho, double dt, double z1, double z2) {
        double dWt = sqrt(dt) * z1;
        double dZt = sqrt(dt) * (rho * z1 + sqrt(1 - rho * rho) * z2); 
    
        double psi = sigma * sqrt(volatilityCurrent) * dWt;
    
        volatilityCurrent += kappa * (theta - volatilityCurrent) * dt + sigma * sqrt(volatilityCurrent) * dZt;
        underlyingCurrent *= exp((riskFreeRate - dividendRate - 0.5 * volatilityCurrent) * dt + sqrt(volatilityCurrent) * (rho * dWt + sqrt(1 - rho * rho) * dZt));
    }

    double lr_americanOption_forGreeks(double St, double K, double T, double r, double sig, double dividend, int N, bool callOption) {
        if (N % 2 == 0) {  // use odd time-steps for correct convergence 
            N += 1;
        }
        double dt = T/N;

        double d1 = d1_(St, K, T, r, sig, dividend);
        double d2 = d2_(St, K, T, r, sig, dividend);

        double pu = peizer_pratt_cdf(d2, N);
        double pd = peizer_pratt_cdf(d1, N);

        double u = exp((r - dividend) * dt) * (pd / pu);
        double d = exp((r - dividend) * dt) * ((1 - pd) / (1 - pu));

        vector<vector<double> > S(N+1, vector<double>(N+1));
        vector<vector<double> > C(N+1, vector<double>(N+1));

        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= i; j++) {
                S[i][j] = St * pow(u, j) * pow(d, i-j);
            }
        }

        for (int j = 0; j <= N; j++) {
            C[N][j] = max(0.0, ((callOption? 1 : -1) * (S[N][j] - K)));
        }

        for (int i = N-1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                C[i][j] = max((exp(-r * dt) * (pu*C[i+1][j+1] + (1-pu)*C[i+1][j])), ((callOption ? 1 : -1) * (S[i][j] - K)));
            }
        }
        return C[0][0];
    }
}

#endif