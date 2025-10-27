# BlackScholes

A fast implementation of the black scholes algorithm in C++

## Build

```bash
cmake -B build/
cmake --build build
```

## Run

```bash
./build/blackscholes
```

## The Black Scholes Model


$$
C_t = N(d_{1})S_t - N(d_2)Ke^{-r\tau}
$$

$$
d_1 = \frac{ln(\frac{S_{t}}{K}) + \tau(r+\frac{\sigma^{2}}{2})}  {\sigma\sqrt{\tau}}
$$

$$
d_{2} = d_{1} - \sigma\sqrt{\tau}
$$

Put price :

$$
\begin{align}
P{t}&=Ke^{-r\tau} - S{t} + C_{t}
\Leftrightarrow P{t}&=N(-d_{2})Ke^{-r\tau} - N(-d_{1})S{t}
\end{align}
$$

* N(x) : The standard normal cumulative distribution function


$$
N(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} e^{-\frac{1}{2}z^{2}}dz
$$



* C(t) : the call option price at time t
* S(t) : the price of the underlying asset at time t
* $\sigma$: the standard deviation of the stock's return (volatility?)
* K : Strike Price (aka exercise price)
* r : annualized the risk-free interest rate
* t : time in years (t_{0} current year)
* T : Time of option expiration
* $\tau$ : the time until maturity ($\tau=(T-t)$)
* N($d_1$) and N($d_2$) are cumulative distribution functions of the standard normal distribution
with

