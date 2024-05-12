#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double binomialEuropeanOption(double S, double K, double T, double r, double u, double d, int N, bool callOption) {
    double dt = T/N;                         // one step change in time
    double p = (exp(r * dt) - d) / (u - d);  // probability of going up
    double q = 1 - p;                        // probability of going down
    double disc = exp(-r * dt);              // discounted 

    double St[N+1]; // asset prices
    double C[N+1];  // option prices

    // initialize the asset prices at maturity (step N) depending on the # of ups (i)
    St[0] = S * pow(d, N);
    for (int i = 1; i <= N; i++) {
        St[i] = St[i - 1] * (u/d);
    }

    // initialize the possible option prices at maturity (step N)
    for (int i = 0; i <= N; i++) {
        C[i] = max(0.0, ((callOption ? 1 : -1) * (St[i] - K)));
    }

    // step back in the tree to get the option price where # of ups is j
    for (int i = N-1; i >= 0; i--) {
        for (int j = 0; j <= i; j++) {
            C[j] = disc * (p * C[j+1] + q * C[j]);
        }
    }

    return C[0];
}

double binomialAmericanOption(double S, double K, double T, double r, double u, double d, int N, bool callOption) {
    double dt = T/N;                         // one step change in time
    double p = (exp(r * dt) - d) / (u - d);  // probability of going up
    double q = 1 - p;                        // probability of going down
    double disc = exp(-r * dt);              // discounted 

    double St[N+1]; // asset prices
    double C[N+1];  // option prices

    // initialize the asset prices at maturity (step N) depending on the # of ups (i)
    St[0] = S * pow(d, N);
    for (int i = 1; i <= N; i++) {
        St[i] = St[i - 1] * (u/d);
    }

    // initialize the possible option prices at maturity (step N)
    for (int i = 0; i <= N; i++) {
        C[i] = max(0.0, ((callOption ? 1 : -1) * (St[i] - K)));
    }

    // step back in the tree to get the option price where # of ups is j
    for (int i = N-1; i >= 0; i--) {
        for (int j = 0; j <= i; j++) {
            C[j] = disc * (p * C[j+1] + q * C[j]);
            St[j] = St[j] / (callOption ? u : d);
            C[j] = max(C[j], ((callOption ? 1 : -1) * (St[j] - K)));  // early exercise test
        }
    }

    return C[0];
}


int main()
{
    double u = 1.2;  // up factor (volatility)
    double d = 1/u;  // down factor (volatility)
    double S = 100;  // current stock price
    double K = 100;  // option strike price
    double r = 0.06;  // risk-free interest rate
    double T = 1;  // time until expiration (in days)
    int N = 3;  // number of time steps
    bool callOption = true;  // true for Call, false for Put
    
    char* optionType = "american";  // write "american" for american options or "european" for european options

    if (optionType == "european") {
        cout << "European " << (callOption ? "Call" : "Put") << " Option Value: " << binomialEuropeanOption(S, K, T, r, u, d, N, callOption) << endl;
    }
    else if (optionType == "american") {
        cout << "American " << (callOption ? "Call" : "Put") << " Option Value: " << binomialAmericanOption(S, K, T, r, u, d, N, callOption) << endl;
    }
    else {
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}