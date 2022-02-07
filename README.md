# A-MFDCCA
codes for A-MFDCCA

-Available in python3

-Please import "numpy" and "math"

-Can implement MF-ADCCA including DFA, DCCA, A-DFA, A-DCCA, MFDFA, MFDCCA, A-MFDFA analysis.

## paper
-S. Kakinaka and K. Umeno, Exploring asymmetric multifractal cross-correlations of price-volatility and asymmetric volatility dynamics in cryptocurrency markets. Physica A 581, 126237 (2021) https://doi.org/10.1016/j.physa.2021.126237

## Usage

If one wants to implement the index-based MF-ADCCA using log-return financial series, then

	basic_dcca(x, y, Q, m=2, trend_base=None, skip_agg=False, observations=100, asymmetry_base='index')

If one wants to implement the return-based A-DFA using financial price series, then

	basic_dcca(x, x, Q=2, m=2, trend_base=None, skip_agg=True, observations=100, asymmetry_base='return')

## Return

S, Fqs, Fhq, S, Fqs_plus, Fhq_plus, S, Fqs_minus, Fhq_minus
- 1-2: estimates for overall trend, 4-5: estimates for positive trend, 7-8: estimates for negative trend
- S: scales, Fqs: fluctuation functions, Fhq: generalized hurst exponent

## Details
	Args:
            x(array(float))  : Target time series data.
            y(array(float))  : Target time series data.
            S(array(int))    : Intervals which divides culmative sum time series. Needs to a positive integer.
            m(int)           : Degree of polynomial fit for each divided segment. Generally, m=2 is recomendend
            Q(array(int))    : fluctuation q-th order.
            trend_base       : designated time series used for defining asymmetric trends.
            skip_agg(bool)   : Whether use cumsum for profile. If series is stationary (ex. log-returns), then "False". If non-stationary (ex. level data), then "True".
                               Its not needed.
            asymmetry_base  : Which criterion used for defining asymmetric trends.
        Returns:
            array(float)     : Return generalized hurst exponent (in np.array) of each Q-th order.
