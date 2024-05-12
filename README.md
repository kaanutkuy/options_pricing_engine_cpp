#### An object-oriented Options Pricing library in C++. Efficiently and accurately price call & put Vanilla and Exotic options and their Greeks.

### Running the Pricer
Build and run the "main.cpp" file, and follow the directions using the command line interface. It will prompt you to create an option by specifying the type of option, and the parameter inputs of the option, and then it will ask for a name to store the option object in memory. After you create it, you can retrieve the option's price, and depending on the type of it, the Greeks as well.

### Types of Options Available
The pricer is able to price Vanilla (European and American) and Exotic (Asian, Barrier, Lookback) options using Black Scholes formula, Binomial/Trinomial trees, and Monte Carlo simulations. Floating Strike is also available for Asian and Lookback exotic options. Pricing of the Greeks is available for European and American options, either with the Black Scholes Greeks modified for yearly dividends for the former, and with finite difference ratios using the Leisen Reimer binomial tree for the latter.
