%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  %%%%%
%%%%      NICTA Energy System Test Case Archive (NESTA) - v0.7.0      %%%%%
%%%%              Optimal Power Flow - Typical Operation              %%%%%
%%%%                         05 - June - 2017                         %%%%%
%%%%                                                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Power flow data for IEEE 13 bus test case.
%  This data was converted from IEEE Common Data Format
%
%  Converted from IEEE CDF file from:
%       http://www.ee.washington.edu/research/pstca/
%
%  CDF Header:
%  IEEE 13 Bus Test Case
%
function mpc = case13_ieee
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	    0.00000	 0.0	 1	    1.05000	    0.95000;
	2	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.04017	   -4.25303	 0.0	 1	    1.05000	    0.95000;
	3	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.00788	  -12.18881	 0.0	 1	    1.05000	    0.95000;
	4	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.01043	   -9.73410	 0.0	 1	    1.05000	    0.95000;
	5	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.01324	   -8.24923	 0.0	 1	    1.05000	    0.95000;
	6	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	  -13.80994	 0.0	 1	    1.05000	    0.95000;
	7	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.04402	  -12.83516	 0.0	 1	    1.05000	    0.95000;
	8	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	  -12.83516	 0.0	 1	    1.05000	    0.95000;
	9	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.04096	  -14.45372	 0.0	 1	    1.05000	    0.95000;
	10	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.03679	  -14.63018	 0.0	 1	    1.05000	    0.95000;
	11	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.04473	  -14.35082	 0.0	 1	    1.05000	    0.95000;
	12	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.04467	  -14.67869	 0.0	 1	    1.05000	    0.95000;
	13	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.03948	  -14.75038	 0.0	 1	    1.05000	    0.95000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	3	 0.0	 0.0	 250	 -250	 1.06	 100.0	 1	 250	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
	5	 0.0	 0.0	 250	 -250	 1.04017	 100.0	 1	 275	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
	7	 0.0	 0.0	 250	 -250	 1.00788	 100.0	 1	 300	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
	11	 0.0	 0.0	 250	 -250	 1.06	 100.0	 1	 200	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
	12	 0.0	 0.0	 250	 -250	 1.06	 100.0	 1	 225	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
	13	 0.0	 0.0	 250	 -250	 1.06	 100.0	 1	 225	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.2	 0.0	 0.0	 0.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   10.0000	   5.000000	   2.00000;
	2	 0.0	 0.0	 3	   20.0000	   10.00000	   4.00000;
	2	 0.0	 0.0	 3	   30.0000	   15.00000	   8.00000;
	2	 0.0	 0.0	 3	   40.0000	   20.00000	   10.0000;
	2	 0.0	 0.0	 3	   50.0000	   25.00000	   5.00000;
	2	 0.0	 0.0	 3	   40.0000	   20.00000	   10.0000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	7	2	0.129469738	0.391477398	0.378788	200	200	200	0	0	1	-30	30	;
	2	3	0.070823886	0.113349152	0.094697	200	200	200	0	0	1	-30	30	;
	2	5	0.12562504	0.128030344	0.094697	200	200	200	0	0	1	-30	30	;
	2	9	0.129469738	0.391477398	0.378788	200	200	200	0	0	1	-30	30	;
	5	6	0.075375024	0.076818206	0.0568182	200	200	200	0	0	1	-30	30	;
	3	4	0.000647349	0.001957387	0.00189394	200	200	200	0	0	1	-30	30	;
	9	13	0.000647349	0.001957387	0.00189394	200	200	200	0	0	1	-30	30	;
	9	11	0.064734869	0.195738699	0.189394	200	200	200	0	0	1	-30	30	;
	9	12	0.075375024	0.076818206	0.0568182	200	200	200	0	0	1	-30	30	;
	13	10	0.075299898	0.0409312	0.094697	200	200	200	0	0	1	-30	30	;
	12	8	0.203409156	0.077636388	0.1515152	200	200	200	0	0	1	-30	30	;
	12	1	0.037761376	0.038281262	0.0568182	200	200	200	0	0	1	-30	30	;
];

% INFO    : === Translation Options ===
% INFO    : Phase Angle Bound:           30.0 (deg.)
% INFO    : Line Capacity Model:         stat
% INFO    : Gen Active Capacity Model:   stat
% INFO    : Gen Reactive Capacity Model: am50ag
% INFO    : Gen Active Cost Model:       stat
% INFO    : AC OPF Solution File:        nesta_case14_ieee.m.opf.sol
% INFO    : Line Capacity PAB:           15.0 (deg.)
% INFO    :
% INFO    : === Generator Classification Notes ===
% INFO    : SYNC   3   -     0.00
% INFO    : NG     2   -   100.00
% INFO    :
% INFO    : === Generator Active Capacity Stat Model Notes ===
% INFO    : Gen at bus 1 - NG	: Pg=232.4, Pmax=332.4 -> Pmax=362   samples: 15
% INFO    : Gen at bus 2 - NG	: Pg=40.0, Pmax=140.0 -> Pmax=63   samples: 1
% INFO    : Gen at bus 3 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    : Gen at bus 6 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    : Gen at bus 8 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    :
% INFO    : === Generator Reactive Capacity Atmost Max 50 Percent Active Model Notes ===
% INFO    : Gen at bus 2 - NG	: Pmax 63.0, Qmin -40.0, Qmax 50.0 -> Qmin -32.0, Qmax 32.0
% INFO    :
% INFO    : === Generator Active Cost Stat Model Notes ===
% INFO    : Updated Generator Cost: NG - 0.0 20.0 0.0430293 -> 0 1.02613295282 0
% INFO    : Updated Generator Cost: NG - 0.0 20.0 0.25 -> 0 0.480811431055 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    :
% INFO    : === Line Capacity Stat Model Notes ===
% WARNING : Missing data for branch flow stat model on line 1-2 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.01938 x=0.05917
% INFO    : Updated Thermal Rating: on line 1-2 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 472
% WARNING : Missing data for branch flow stat model on line 1-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05403 x=0.22304
% INFO    : Updated Thermal Rating: on line 1-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 128
% WARNING : Missing data for branch flow stat model on line 2-3 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.04699 x=0.19797
% INFO    : Updated Thermal Rating: on line 2-3 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 145
% WARNING : Missing data for branch flow stat model on line 2-4 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05811 x=0.17632
% INFO    : Updated Thermal Rating: on line 2-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 158
% WARNING : Missing data for branch flow stat model on line 2-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05695 x=0.17388
% INFO    : Updated Thermal Rating: on line 2-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 161
% WARNING : Missing data for branch flow stat model on line 3-4 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.06701 x=0.17103
% INFO    : Updated Thermal Rating: on line 3-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 160
% WARNING : Missing data for branch flow stat model on line 4-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.01335 x=0.04211
% INFO    : Updated Thermal Rating: on line 4-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 664
% WARNING : Missing data for branch flow stat model on line 4-7 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.20912
% INFO    : Updated Thermal Rating: on transformer 4-7 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 141
% WARNING : Missing data for branch flow stat model on line 4-9 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.55618
% INFO    : Updated Thermal Rating: on transformer 4-9 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 53
% WARNING : Missing data for branch flow stat model on line 5-6 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.25202
% INFO    : Updated Thermal Rating: on transformer 5-6 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 117
% WARNING : Missing data for branch flow stat model on line 6-11 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.09498 x=0.1989
% INFO    : Updated Thermal Rating: on line 6-11 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 134
% WARNING : Missing data for branch flow stat model on line 6-12 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.12291 x=0.25581
% INFO    : Updated Thermal Rating: on line 6-12 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 104
% WARNING : Missing data for branch flow stat model on line 6-13 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.06615 x=0.13027
% INFO    : Updated Thermal Rating: on line 6-13 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 201
% WARNING : Missing data for branch flow stat model on line 7-8 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.17615
% INFO    : Updated Thermal Rating: on line 7-8 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 167
% WARNING : Missing data for branch flow stat model on line 7-9 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.11001
% INFO    : Updated Thermal Rating: on line 7-9 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 267
% WARNING : Missing data for branch flow stat model on line 9-10 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.03181 x=0.0845
% INFO    : Updated Thermal Rating: on line 9-10 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 325
% WARNING : Missing data for branch flow stat model on line 9-14 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.12711 x=0.27038
% INFO    : Updated Thermal Rating: on line 9-14 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 99
% WARNING : Missing data for branch flow stat model on line 10-11 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.08205 x=0.19207
% INFO    : Updated Thermal Rating: on line 10-11 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 141
% WARNING : Missing data for branch flow stat model on line 12-13 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.22092 x=0.19988
% INFO    : Updated Thermal Rating: on line 12-13 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 99
% WARNING : Missing data for branch flow stat model on line 13-14 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.17093 x=0.34802
% INFO    : Updated Thermal Rating: on line 13-14 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 76
% INFO    :
% INFO    : === Voltage Setpoint Replacement Notes ===
% INFO    : Bus 1	: V=1.06, theta=0.0 -> V=1.06, theta=0.0
% INFO    : Bus 2	: V=1.045, theta=-4.98 -> V=1.04017, theta=-4.25303
% INFO    : Bus 3	: V=1.01, theta=-12.72 -> V=1.00788, theta=-12.18881
% INFO    : Bus 4	: V=1.019, theta=-10.33 -> V=1.01043, theta=-9.7341
% INFO    : Bus 5	: V=1.02, theta=-8.78 -> V=1.01324, theta=-8.24923
% INFO    : Bus 6	: V=1.07, theta=-14.22 -> V=1.06, theta=-13.80994
% INFO    : Bus 7	: V=1.062, theta=-13.37 -> V=1.04402, theta=-12.83516
% INFO    : Bus 8	: V=1.09, theta=-13.36 -> V=1.06, theta=-12.83516
% INFO    : Bus 9	: V=1.056, theta=-14.94 -> V=1.04096, theta=-14.45372
% INFO    : Bus 10	: V=1.051, theta=-15.1 -> V=1.03679, theta=-14.63018
% INFO    : Bus 11	: V=1.057, theta=-14.79 -> V=1.04473, theta=-14.35082
% INFO    : Bus 12	: V=1.055, theta=-15.07 -> V=1.04467, theta=-14.67869
% INFO    : Bus 13	: V=1.05, theta=-15.16 -> V=1.03948, theta=-14.75038
% INFO    : Bus 14	: V=1.036, theta=-16.04 -> V=1.02209, theta=-15.60843
% INFO    :
% INFO    : === Generator Setpoint Replacement Notes ===
% INFO    : Gen at bus 1	: Pg=232.4, Qg=-16.9 -> Pg=208.317, Qg=0.0
% INFO    : Gen at bus 1	: Vg=1.06 -> Vg=1.06
% INFO    : Gen at bus 2	: Pg=40.0, Qg=42.4 -> Pg=63.0, Qg=26.964
% INFO    : Gen at bus 2	: Vg=1.045 -> Vg=1.04017
% INFO    : Gen at bus 3	: Pg=0.0, Qg=23.4 -> Pg=0.0, Qg=29.63
% INFO    : Gen at bus 3	: Vg=1.01 -> Vg=1.00788
% INFO    : Gen at bus 6	: Pg=0.0, Qg=12.2 -> Pg=0.0, Qg=13.23
% INFO    : Gen at bus 6	: Vg=1.07 -> Vg=1.06
% INFO    : Gen at bus 8	: Pg=0.0, Qg=17.4 -> Pg=0.0, Qg=9.618
% INFO    : Gen at bus 8	: Vg=1.09 -> Vg=1.06
% INFO    :
% INFO    : === Writing Matpower Case File Notes ===
