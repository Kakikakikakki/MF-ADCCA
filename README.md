# MF-ADCCA
codes for MF-ADCCA: MF-ADCCA_code.py

- Available in python3

- Please import `numpy` and `math`

- Can implement MF-ADCCA including DFA, DCCA, A-DFA, A-DCCA, MFDFA, MFDCCA, A-MFDFA analysis.

## paper
when using this code, please cite the following paper:

> S. Kakinaka and K. Umeno, Exploring asymmetric multifractal cross-correlations of price-volatility and asymmetric volatility dynamics in cryptocurrency markets. Physica A 581, 126237 (2021) https://doi.org/10.1016/j.physa.2021.126237

## Usage
Note that in asymmetric analysis, the criterion of which series you use to define up and down trends may be cruical.

The `asymmetry_base` can be set as :
- `index`, when you want the rise and fall of index itself to be the criterion to define up and down trends.
- `return`, when you want the rise and fall of return series to be the criterion to define up and down trends.
- `optional`, when you want to use an alternative series to be the criterion to define up and down trends.

There are three patterns of asymmetric versions related to DFA-based analysis:

If one wants to implement the index-based MF-ADCCA using financial return series (log-returns), then

	basic_dcca(x, y, Q=2, m=2, trend_base=None, skip_agg=False, observations=100, asymmetry_base='return')

If one wants to implement the return-based A-DFA, substitue y to x:

	basic_dcca(x, x, Q=2, m=2, trend_base=None, skip_agg=True, observations=100, asymmetry_base='return')

If one wants to designate a specific series for the benchmark aganist defining asymmetry, then by using your own series at hand (my_series),
	
	basic_dcca(x, x, Q=2, m=2, trend_base=my_series, skip_agg=True, observations=100, asymmetry_base='optional')

If one is intrested only in symmetric analysis (for example, MFDCCA), benefit from the MF-ADCCA estimates, that is, utilize the estimates of the overall trend:
	
	basic_dcca(x, y, Q=2, m=2, trend_base=my_series, skip_agg=True, observations=100, asymmetry_base='return')[0:2]

## Return

`list: [S, Fqs, Fhq, S, Fqs_plus, Fhq_plus, S, Fqs_minus, Fhq_minus]`
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
            asymmetry_base   : Which criterion used for defining asymmetric trends.
    Returns:
            array(float)     : list of estimates e.g., Return generalized hurst exponent (in np.array) of each Q-th order.

## Further parameter settings for estimation

You can rearrange the parameters necessary for implmenting the analysis.
- observations: the number of observational scale points used for the log-log linear fit when estimating (genralized) Hurst exponents (default set at 100).
- scales are set at: `s_min = max(20, int(np.floor(N/100)))` and `s_max = min(20*s_min, int(np.floor(N/10)))`
  
