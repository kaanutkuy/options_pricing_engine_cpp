#### An object-oriented Options Pricing library in C++. Efficiently and accurately price call & put Vanilla and Exotic options and their Greeks.

### Running the Pricer
Build and run the "main.cpp" file, and follow the directions using the command line interface. It will prompt you to create an option by specifying the type of option, and the parameter inputs of the option, and then it will ask for a name to store the option object in memory. After you create it, you can retrieve the option's price, and depending on the type of it, the Greeks as well.

### Types of Options Available
The pricer is able to price Vanilla (European and American) and Exotic (Asian, Barrier, Lookback) options using Black Scholes formula, Binomial/Trinomial trees, and Monte Carlo simulations. Pricing Floating Strike options is also available for Asian and Lookback exotic options. Pricing of the Greeks is available for European and American options, either with the Black Scholes Greeks modified for yearly dividends for the former, and with finite difference ratios using the Leisen Reimer binomial tree nodes for the latter.

### Project Structure
| File Name | Description |
| ---- | ---- |
| option_class.hpp | Option Class template. Contains "Option" and "Floating Option" parent classes to create, modify, and display various vanilla & exotic options with fixed or floating strikes. Includes getters for input parameters of the underlying, Greeks and the option price, setters for input parameters of the underlying, and a display function used to display the option's qualities in the main file. |
| utility_functions.hpp | Contains all the neccessary functions like distribution functions, Greeks calculations helpers, forward prices, etc. that the pricers use in the pipeline. Includes probability and cumulative density functions, forward price functions, Black-Scholes Greek delta, gamma, vega, theta, rho functions modified for yearly dividend return, sign function, Peizer-Pratt cumulative distribution function for calculating up & down probabilities in the Leisen-Reimer binomial-tree American option pricer, a Heston stochastic volatility model helper function that updates the next underlying price & volatility, and a Leiser-Reimer binomial-tree pricer to calculate elements of the finite differences formula in American options Greeks.|
| vanilla_options.cpp | Holds the European and American option classes that inherit from the parent Option class. Includes functions to price the options and its Greeks, specifically, a Black-Scholes pricer modified for yearly dividend returns, Monte Carlo simulation pricers for constant and stochastic volatility that use antithetic and delta,gamma control variates, Greeks calculation function for European options, and Leisen-Reimer binomial-tree pricer, Cox-Ross-Rubinstein trinomial-tree pricer, and Greeks calculation function for American options.|
| exotic_options.cpp | Holds the Asian, Lookback and Barrier exotic options that inherit from parent Option and Floating-strike option classes. Includes functions to price the options, all using Monte Carlo simulations. |
| main.cpp | Main file to run the program on. Contains functions for a user to interact with the command-line-interface, create and price available options. Can be modified for different capabilities by calling the various functions of the pricers. | 

## Pricers' Functionalities
### Black-Scholes European Option Pricer:
The _blackScholes_europeanOption_ function in the EuropeanOption class in vanilla_options.cpp calculates the price of a European option (call or put) using the Black-Scholes formula, incorporating dividends. The function calculates the price of the option by:

- For **Call** Options:
  C = $e^{-rT} \left( S e^{(r - q)T} N(d1) - K N(d2) \right)$
- For **Put** Options:
  P = $e^{-rT} \left( K N(-d2) - S e^{(r - q)T} N(-d1) \right)$

where:
- S: Underlying price
- K: Strike price
- T: Time to expiry in years
- r: Risk-free rate
- σ: Volatility
- q: Yearly Dividend rate

and _N(x)_ is the cumulative distribution function of the standard normal distribution that represents the probability the option will expire in-the-money, d1 = $\frac{\ln\left(\frac{S}{K}\right) + (r - q + 0.5\sigma^2)T}{\sigma \sqrt{T}}$ and d2 = $\frac{\ln\left(\frac{S}{K}\right) + (r - q - 0.5\sigma^2)T}{\sigma \sqrt{T}}$ represent the "moneyness" of the option where d1 is how many standard deviations the log of the price ratio $S/K$ is from the option's payoff at expiration adjusted for the risk-free rate and dividends and d2 adjusts d1 down to reflect the drift rate of the stock under the risk-neutral measure.

### Monte Carlo Simulations for Exotic & European Options:
The library computes the prices of the Exotic options (Asian, Lookback, and Barrier) using Monte Carlo simulations as these options are path dependent and have complex payoff structures that depend on the underlying's path during its lifetime. The library also gives the user the ability to choose to price European options with Monte Carlo if they wanted to, as these simulations can represent thousands of "possible" movements of the underlying's price and parameters that might represent the real-world conditions better depending on the behaviour of the input parameters.

The Monte Carlo simulations in this library use two methods to **reduce variance**: Antithetic variates, and Delta & Gamma control variates. By using these two known methods, the pricers require less simulations to regress to the actual price accurately, saving total computation time.

Key parameters of the simulations are:
- St: Underlying asset price at time t
- K: Strike price
- T: Time to expiry
- r: Risk-free rate
- q: Yearly Dividend rate
- σ: Volatility
- Δt: Time step $\Delta t = \frac{T}{\text{timeSteps}}$
- ν: Drift term $\nu = r - q - 0.5 \sigma^2$
- Exponential growth factor adjusted for dividends: $e^{r \Delta t - q \Delta t}$
- Auxiliary term for gamma control variate $e_{\gamma} = e^{(2(r - q) + \sigma^2)\Delta t} - 2e^{r \Delta t - q \Delta t} + 1$
- Beta coefficients for delta and gamma control variates: $\beta_1$, $\beta_2$ (came from precomuted regression)

#### Simulation Process:
1. **Initialization**
  - initializes the constants & equations
  - uses Mersenne Twister engine (a pseudo-random normal number generator library) to create a normal distribution to generate numbers for the underlying's increments for every timestep
    - Note: might be updated to a quasi-random number generator (like Sobol) to create a more homogenous, less clustered distribution 
2. **Simulation Loops**
  - For each simulation, antithetic and control variates are set to initial prices
    - For each timestep of the simulation:
      
      
    

3. **Option Payoffs**
  - 

### Leisen-Reimer Binomial Tree American Option Pricer:

### Cox-Ross-Rubinstein Trinomial Tree American Option Pricer:



