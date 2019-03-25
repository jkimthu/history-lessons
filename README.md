# history-lessons
Probing cell division under changing environments


A. Mean interdivision time vs. mean growth rate
   No deviation from expected curve in fluc conditions
   - note: currently interdiv time is cleaned by removing outliers beyond 3 std
   - note: growth rate is not cleaned up in a similar way, however

A2. Mean birth size vs. mean growth rate
   No deviation from expected curve in fluc conditions
   - note: currently interdiv time is cleaned by removing outliers beyond 3 std
   - note: growth rate is not cleaned up in a similar way, however


B1. Mean division size vs mean birth size, population level


B2. Division size vs. birth size, single-cell (heated scatter plot)
	- trimmed for steady-state (post 3 hrs) based on birth time
	- for prior to steady-state, trimming by birth time gives 2x
	- for prior to steady-state, trimming by division time (< 2.5h) gives....


C1. Traits associated with being above or below average in interdivision size
	- isolates two populations from each condition: (1) one std larger than average division size and (2) one std smaller than average division size
	- compares traits between the two subpopulations for all conditions

	- main findings:
	i. interdivions are longer in cells that divide larger (fluc), whereas in steady interdivs ~ equal
	ii. mean growth rate is higher in cells that divide smaller (fluc), whereas in steady ~ equal
	iii. both of these are consistent with steady-state relationships between growth rate and interdivision time
	-- overall: this represents a regulation at the single cell level


C2. Traits binned by nutrient phase at birth
	- nutrient phase defines the start of the low as the first bin (bin 1) and the start of the high being the the start of the second half of bins (if 24 bins total, high phase starts at bin 13) 



D. Bootstrap hypothesis testing of test statistics between a single lineage and random subsamples from the whole population


E & F. assess traits by nutrient signal classifications
		- classifications are performed by sectioning signal into 5ths
		- these 5ths are assigned "high" or "low" based on the majority signal
		- the binary sequence of these fifths are given a class number, by which cell cycles are binned
		- the cell cycle stats or "traits" are measured then for each class bin

		Program E does one condition per designated set of experiments and produces a plot for each experiment, of trait mean and stdev by class.

		Program F compiles cell cycles from each listed experiment and produces a single plot for all data. This plot compares signals by splitting classifications into two more data points:
			1. fraction of signal spent in high
			2. signal type: low to high, high to low, high-low-high, low-high-low

		Both programs, given the manner in which classifications are defined, only function with data from 60 min periods. A new classification scheme needs to be developed to consider other data, including 15 min periods.



G. size vs time (at birth, at division, added)
	- what is the timescales at which cell size stabilizes? how does it compare with that of mu?
	- to get at this, four plots:
	1. birth size vs time at birth
	2. division size vs time at division
	3. added size vs time at birth
	4. added size vs time at division


