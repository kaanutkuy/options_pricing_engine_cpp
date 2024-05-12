#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>

#include "utility_functions.hpp"

using namespace std;
using namespace std::chrono;


// Heston model dynamics
void hestonDynamics(double& St, double& vt, double r, double dividend, double kappa, double theta, double sigma, double rho, double dt, double z1, double z2) {
    double dWt = sqrt(dt) * z1;
    double dZt = sqrt(dt) * (rho * z1 + sqrt(1 - rho * rho) * z2); 
    
    double psi = sigma * sqrt(vt) * dWt;
    
    vt += kappa * (theta - vt) * dt + sigma * sqrt(vt) * dZt;
    St *= exp((r - dividend - 0.5 * vt) * dt + sqrt(vt) * (rho * dWt + sqrt(1 - rho * rho) * dZt));
}


// American options pricing
double cox_ross_rubinstein_trinomial(double St, double K, double T, double r, double sig, double dividend, int N, bool callOption) {
    double dt = T/N;

    double nu = (r - dividend - 0.5*sig*sig);
    double lambda = sqrt(3/2);  // Kamrad and Ritchken best convergence rate sqrt(3/2)
    double pu = 1/(2*lambda*lambda) + 0.5 * (nu / (lambda*sig)) * sqrt(dt);
    double pd = 1/(2*lambda*lambda) - 0.5 * (nu / (lambda*sig)) * sqrt(dt);
    double pm = 1 - pu - pd;

    double u = exp(lambda * sig * sqrt(dt));
    double p = 1/u;

    double S[N+1][2*N+1];
    double C[N+1][2*N+1];
    // vector<vector<double> > S(N+1, vector<double>(2*N+1));
    // vector<vector<double> > C(N+1, vector<double>(2*N+1));

    for (int i = N; i >= 0; i--) {
        for (int j = -i; j <= i; j++) {
            S[i][j] = St * pow(u, j);
        }
    }

    for (int j = N; j >= -N; j--) {
        C[N][j] = max(0.0, ((callOption ? 1 : -1) * (S[N][j] - K)));
    }

    for (int i = N-1; i >= 0; i--) {
        for (int j = i; j >= -i; j--) {
            C[i][j] = max((exp(-r*dt) * (pu*C[i+1][j+1] + pm*C[i+1][j] + pd*C[i+1][j-1])), ((callOption ? 1 : -1) * (S[i][j] - K)));
        }
    }
    //double delta = (C[1][1] - C[1][0]) / (S[1][1] - S[1][0]);
    //cout << C[1][1] << endl;
    //cout << C[1][0] << endl;
    //cout << S[1][1] << endl;
    //cout << S[1][0] << endl;
    //cout << delta << endl;
    return C[0][0];
}

// Leisen-Reimer model for pricing American options
double leisen_reimer_binomial(double St, double K, double T, double r, double sig, double dividend, int N, bool callOption) {
    if (N % 2 == 0) {  // use odd time-steps for correct convergence 
        N += 1;
    }
    double dt = T/N;

    double d1 = UtilityFunctions::d1_(St, K, T, r, sig, dividend);
    double d2 = UtilityFunctions::d2_(St, K, T, r, sig, dividend);

    double pu = UtilityFunctions::peizer_pratt_cdf(d2, N);
    double pd = UtilityFunctions::peizer_pratt_cdf(d1, N);

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

    double delta = (C[1][1] - C[1][0]) / (S[1][1] - S[1][0]);
    double gamma = (((C[2][2] - C[2][1]) / (S[2][2] - S[2][1])) - ((C[2][1] - C[2][0]) / (S[2][1] - S[2][0]))) / (0.5 * (S[2][2] - S[2][0]));

    cout << "delta " << delta << endl;
    cout << "gamma " << gamma << endl;

    return C[0][0];
}

//double leisen_reimer_heston_SV(double St, double K, double T, double r, double sig, double dividend, int N, bool callOption) {
//    if (N % 2 == 0) {  // use odd time-steps for correct convergence 
//        N += 1;
//    }
//    double dt = T/N;
//
//    double d1 = d1_(St, K, T, r, sig, dividend);
//    double d2 = d2_(St, K, T, r, sig, dividend);
//
//    double pu = peizer_pratt_cdf(d2, N);
//    double pd = peizer_pratt_cdf(d1, N);
//
//    double u = exp((r - dividend) * dt) * (pd / pu);
//    double d = exp((r - dividend) * dt) * ((1 - pd) / (1 - pu));
//
//    vector<vector<double> > S(N+1, vector<double>(N+1));
//    vector<vector<double> > C(N+1, vector<double>(N+1));
//
//    for (int i = 0; i <= N; i++) {
//        for (int j = 0; j <= i; j++) {
//            S[i][j] = St * pow(u, j) * pow(d, i-j);
//        }
//    }
//
//    for (int j = 0; j <= N; j++) {
//        C[N][j] = max(0.0, ((callOption ? 1 : -1) * (S[N][j] - K)));
//    }
//
//    for (int i = N-1; i >= 0; i--) {
//        for (int j = 0; j <= i; j++) {
//            C[i][j] = max((exp(-r*dt) * (pu*C[i+1][j+1] + (1 - pu)*C[i+1][j])), ((callOption ? 1 : -1) * (S[i][j] - K)));
//
//            // recombine option prices to improve computational efficiency (to not have exploding number of prica calculations caused by changing volatility)
//            if (j > 0) {
//                double C_down = exp(-r * dt) * (pu*C[i + 1][j] + (1 - pu)*C[i + 1][j - 1]);  // recombine up and down option prices
//                C[i][j] = max(C[i][j], C_down);
//            }
//
//        }
//    }
//    return C[0][0];
//}


int main() {
    double sig = 0.15;  // volatility 
    double dividend = 0;  // dividend pay per share
    double S = 33.75;  // current stock price
    double K = 35;  // option strike price
    double r = 0.055;  // risk-free interest rate
    double T = 0.75;  // time until expiration (in years)
    int N = 10;  // number of time steps (>= 100 steps for best performance)
    bool callOption = false;  // true for Call, false for Put

    auto start_crr = high_resolution_clock::now();
    //cout << "CRR American " << (callOption ? "Call" : "Put") << " Option Value: " << cox_ross_rubinstein_trinomial(S, K, T, r, sig, dividend, N, callOption) << endl;
    auto stop_crr = high_resolution_clock::now();
    auto duration_crr = duration_cast<microseconds>(stop_crr - start_crr);
    
    auto start_lr = high_resolution_clock::now();
    //cout << "LR American " << (callOption ? "Call" : "Put") << " Option Value: " << leisen_reimer_binomial(S, K, T, r, sig, dividend, N/2, callOption) << endl;
    auto stop_lr = high_resolution_clock::now();
    auto duration_lr = duration_cast<microseconds>(stop_lr - start_lr);
    cout << leisen_reimer_binomial(S, K, T, r, sig, dividend, N/2, callOption) << endl;
    //cout << "CRR compute time: " << duration_crr.count() << " microseconds" << endl;
    //cout << "LR compute time: " << duration_lr.count() << " microseconds" << endl;

    double dsig = sig * 0.0001;
    double dr = r * 0.0001;
    double dT = T * 0.0001;
    double C1 = leisen_reimer_binomial(S, K, T, r, sig + dsig, dividend, N/2, callOption);
    double C2 = leisen_reimer_binomial(S, K, T, r, sig - dsig, dividend, N/2, callOption);
    double vega = (C1 - C2) / (2*dsig);
    C1 = leisen_reimer_binomial(S, K, T, r + dr, sig, dividend, N/2, callOption);
    C2 = leisen_reimer_binomial(S, K, T, r - dr, sig, dividend, N/2, callOption);
    double rho = (C1 - C2) / (2*dr);
    C1 = leisen_reimer_binomial(S, K, T + dT, r, sig, dividend, N/2, callOption);
    C2 = leisen_reimer_binomial(S, K, T - dT, r, sig, dividend, N/2, callOption);
    double theta = (C1 - C2) / (2*dT);

    //cout << "delta " << delta << endl;
    //cout << "gamma " << gamma << endl;
    cout << "vega " << vega << endl;
    cout << "rho " << rho << endl;
    cout << "theta " << theta << endl;
    

    return EXIT_SUCCESS;
}
