### An object-oriented Options Pricing library in C++. Efficiently and accurately price call & put Vanilla and Exotic options and their Greeks with advanced Monte Carlo Simulations, Binomial/Trinomial Trees and closed-form Black-Scholes.

## Running the Pricer
Build and run the "main.cpp" file, and follow the directions using the command line interface. It will prompt you to create an option by specifying the type of option, and the parameter inputs of the option, and then it will ask for a name to store the option object in memory. After you create it, you can retrieve the option's price, and depending on the type of it, the Greeks as well.

## Types of Options Available
The pricer is able to price Vanilla (European and American) and Exotic (Asian, Barrier, Lookback) options using Black Scholes formula, Binomial/Trinomial trees, and Monte Carlo simulations. Pricing Floating Strike options is also available for Asian and Lookback exotic options. Pricing of the Greeks is available for European and American options, either with the Black Scholes Greeks modified for yearly dividends for the former, and with finite difference ratios using the Leisen Reimer binomial tree nodes for the latter.

## Project Structure
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
The library computes the prices of the Exotic options (Asian, Lookback, and Barrier) using Monte Carlo simulation functions as these options are path dependent and have complex payoff structures that depend on the underlying's path during its lifetime. The library also gives the user the ability to choose to price European options with Monte Carlo if they wanted to, as these simulations can represent thousands of "possible" movements of the underlying's price and parameters that might represent the real-world conditions better depending on the behaviour of the input parameters.

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
- N: number of simulations

#### Simulation Process:
1. **Initialization**
  - initializes the constants & equations
  - uses Mersenne Twister engine (a pseudo-random normal number generator library) to create a normal distribution to generate numbers for the underlying's increments for every timestep
    - Note: might be updated to a quasi-random number generator (like Sobol) to create a more homogenous, less clustered distribution 
2. **Simulation Loops**
  - For each simulation, antithetic and control variates are set to initial prices.
    - For each timestep of the simulation:
      - Generates a random number from a standard normal distribution
      - Calculates the delta and gamma Greeks for the underlying and its negatively correlated antithetic variate
      - Updates the underlying price and its antithetic with: $$S_{t+1} = S_t \exp((r-q-\frac{1}{2} \sigma^2)\Delta t + \sigma \sqrt{\Delta t} Z_{t+1})$$ where $Z_{t+1}$ is the number from the normal distribution
      - Calculates the delta & gamma control variates
    - Takes the average of the antithetic and delta & gamma aided option payoffs to estimate the final option payoff for each simulation
      - **NOTE**: The option payoff calculations differ in the exotic options. For example, Barrier options add to the option sums depending on if a certain threshold is crossed or not, or Lookback options have an inherent exercising function where the option payoff depends on the maximum or minimum value the underlying has been through its lifetime, or Asian options add to the option sums the difference of the final underlying price and its arithmetic/geometric average of this underlying through its lifetime.

3. **Option Payoffs**
  - Aggregates the option payoffs across all simulations, take their average, and discount it back to present value to estimate the option price with: $$\text{Option Price} = e^{-rT} \cdot \frac{1}{N} \cdot \sum_{n=1}^{N} {Payoff_n}$$
  
### Leisen-Reimer Binomial Tree American Option Pricer:
The main pricer of the American option the library uses is the _ls_americanOption_ Leisen Reimer binomial tree function in the AmericanOption class of vanilla_options.cpp file. The reason that it is the preferred pricer of American options is that Leisen Reimer trees require much less timesteps to accurately compute American options compared to other methods, which equates to faster computing of the option price. Also, Leisen Reimer trees assume that the binomial tree is centered around the option's strike price at its expiration, instead of its underlying price, which in turn makes it more suitable if the user believes that the option's strike price should be the main driver of its price instead of the current underlying.

The key parameters of the tree are:
- S: Underlying price
- K: Strike price
- T: Time to expiry in years
- r: Risk-free rate
- σ: Volatility
- q: Yearly Dividend rate
- timesteps: The timesteps to discretize the underlying's and option's lifetime as a grid

#### Tree Structure

1. Leisen Reimer tree first makes sure that the timesteps of the option lifetime is **odd**, because an odd number of timesteps ensures that the final time step, which corresponds to the expiration of the option, is even, helping to maintain symmetry in the tree structure.
2. Calculates d1 & d2 (the moneyness of the option)
3. Calculates up (pu) and down (pd) probabilities that the option can take for the next step of the grid, using the _Peizer-Pratt inversion_ cumulative distribution function.
4. Calculates the magnitude of the up and down movements using the rate of the probabilities of these movements over each other.
5. Calculates the possible values the underlying can take depending on the up & down moves
6. Calculates the final payoffs of the option at the end of its life.
7. Calculates the payoffs of the previous steps by stepping back through the tree grid from the final payoffs, and finds the initial value of the payoff
8. Calculates the delta & gamma Greeks using finite difference ratios from the up & down option and underlying prices of the next steps in the binomial tree

## Current Limitations of the Library & Possible Future Additions

