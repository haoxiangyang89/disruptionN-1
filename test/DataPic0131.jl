#N,K,CT,DT, CTB, DTB,	V_lb 	V_ub
#Cont_Tech_ind values has value 0 and >=1 for each type of technology
#Node# NodeName   Slots   Cont_Tech_type_ind   Disc_Tech_Ind   Cont_Tech_Bat   Disc_Tech_Bat   Volt_lb   Volt_ub
busmat =[
     1     611	    1	            1	              0	              1	              0	              0.95     1.05;
     2     632	    1	            0	              0	              0	              0	              0.95     1.05;
     3     633	    1	            0	              0	              0	              0	              0.95     1.05;
     4     634	    1	            0	              0	              0	              0	              0.95     1.05;
     5     645	    3	            0	              1	              0	              0	              0.95     1.05;
     6     646	    1	            0	              0	              0	              0	              0.95     1.05;
     7     650	    1	            0	              1	              0	              0	              0.95     1.05;
     8     652	    1	            1	              1	              0	              0	              0.95     1.05;
     9     671	    1	            0	              0	              0	              0	              0.95     1.05;
     10    675	    1	            1	              0	              0	              0	              0.95     1.05;
     11    680	    1	            0	              0	              0	              0	              0.95     1.05;
     12    684	    1	            0	              0	              0	              0	              0.95     1.05;
     13    692	    1	            0	              0	              0	              0	              0.95     1.05;
]
#Time periods
 timeperiod =collect(1:numtimeperiods)'

#Active Demand
##Node# NodeName  T1 T2
  Pit =[
  1  611 0.021	0	0	0.03	0.081	0	0	0.028	0	0.088	0	0	0	0	0	0	0.091	0.029	0	1.757	0	0	0.003	0.11	0.11	0	0	0	0	0	0.093	0.03	0	0	0	0	0	0.052	0.074	0	0	0	0	0	0	0.075	0.057	0	0.152	1.597	0	0	0	0.086	0.045	0	0	0	0	0	0.01	0.111	0.013	0	0	0	0	0.04	0.094	0	0	0	0	0	0	0	0.048	0	1.687	0	0	0	0.04	0.088	0	0	0	0	0	0.097	0.035	0	0	0	0	0;
  2  632 0.24	1.432	2.196	0.266	0.266	0.199	0.204	0.197	0.206	0.201	0.202	0.197	0.254	0.202	0.198	0.196	0.197	0.106	1.184	0.135	0.105	0.106	0.106	0.16	0.16	0.12	0.115	0.113	0.036	0.019	0.02	0.084	0.038	0.038	0.02	0.02	0.125	0.108	0.106	0.193	0.196	0.538	0.363	0.324	0.27	0.649	0.435	0.266	0.271	0.27	0.308	0.375	0.369	0.368	0.34	0.268	0.382	0.362	0.332	1.152	0.806	0.268	0.351	2.255	4.097	2.698	2.698	4.332	3.27	4.274	2.052	2.052	2.32	3.056	3.056	1.796	2.33	3.17	2.18	2.89	1.71	1.856	1.788	2.384	0.542	1.84	0.506	0.51	1.681	0.576	0.516	1.978	1.1	1.387	1.387	1.571;
  3  633 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  4  634 1.46	1.307	3.055	1.353	2.168	1.256	1.323	1.332	1.218	1.318	1.184	1.186	1.25	1.394	1.338	2.482	1.37	1.28	1.382	1.92	2.045	1.988	1.988	1.881	1.891	1.952	2.015	1.924	1.878	2.74	1.978	1.978	1.892	1.875	2.305	3.349	4.274	3.321	6.764	7.126	7.752	7.77	7.77	7.86	3.496	3.041	1.98	1.954	1.854	1.852	1.934	2.986	1.926	1.86	2.506	2.256	1.879	5.079	5.392	2.57	2.994	3.364	4.428	2.152	2.18	3.078	3.078	3.376	3.392	3.462	3.434	3.434	3.366	3.392	3.392	3.461	4.472	4.34	4.154	2.507	3.456	3.445	4.274	3.52	3.514	3.406	3.31	3.308	3.3	4.81	3.201	3.266	3.274	3.18	3.18	3.336;
  5  645 0.631	0.614	0.556	0.56	0.554	0.242	0.554	0.365	0.464	0.352	0.423	0.352	0.274	0.232	0.228	0.23	0.226	0.232	0.332	0.349	0.281	0.23	0.228	0.228	0.228	0.226	0.23	0.222	0.33	0.346	0.28	0.226	0.236	0.236	0.229	0.227	0.226	0.224	0.227	0.304	0.394	0.364	0.26	0.24	0.24	0.215	0.21	0.216	0.332	0.334	0.236	0.217	0.193	0.186	0.182	0.188	0.221	0.427	0.568	0.521	0.375	0.656	0.662	0.691	0.432	0.285	0.256	0.366	0.366	0.36	0.356	0.23	0.23	1.209	1.214	1.214	0.486	0.176	0.174	0.176	0.172	0.174	0.377	0.43	0.418	0.316	0.299	0.297	1.092	1.292	1.174	1.042	1.438	1.088	0.954	0.954;
  6  646 0.039	0.037	0.08	0.139	0.131	0.038	0.04	0.039	0.036	0.038	0.04	0.038	0.038	0.048	0.111	0.137	0.096	0.038	0.035	0.038	0.04	0.036	0.04	0.039	0.039	0.04	0.052	0.138	0.135	0.051	0.04	0.036	0.04	0.04	0.04	0.036	0.038	0.038	0.064	0.42	0.64	1.356	1.252	1.165	1.183	1.163	1.125	1.028	0.832	0.768	0.817	0.663	0.122	0.051	0.098	0.197	0.112	0.027	0.057	0.027	0.085	0.13	0.122	0.056	0.028	0.028	0.024	0.024	0.037	0.131	0.124	0.126	0.126	0.026	0.024	0.024	0.03	0.025	0.125	0.122	0.091	0.025	0.024	0.028	0.024	0.026	0.023	0.101	0.124	0.673	1.143	1.233	1.016	1.328	1.328	0.335;
  7  650 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  8  652 0.48	0.379	0.469	2.685	4.649	0.544	3.178	0.679	1.96	0.643	0.907	0.685	0.524	0.787	0.523	0.694	0.617	0.684	0.377	0.747	0.686	0.507	0.749	0.549	0.549	0.726	0.524	0.682	0.571	0.722	0.511	0.564	0.726	0.726	0.642	0.654	0.503	0.682	0.602	0.292	0.471	0.168	0.12	0.138	0.109	0.134	0.134	0.228	0.168	0.157	0.064	0.096	0.106	0.106	0.047	0.047	0.144	0.069	0.046	0.046	0.071	0.081	0.047	0.088	0.046	0.12	0.12	0.177	0.119	0.119	0.133	0.072	0.048	0.15	0.051	0.048	0.155	0.047	0.048	0.156	0.047	0.154	0.045	0.048	0.047	0.138	0.046	0.126	0.073	0.045	0.045	0.045	0.047	0.143	0.046	0.142;
  9  671 3.002	3.002	1.094	0.932	0.621	0.82	0.539	0.799	0.438	0.587	0.646	0.683	0.646	1.23	3.182	2.438	2.94	1.163	1.769	0.708	0.687	0.687	0.556	0.308	0.299	0.337	0.42	0.372	0.448	1.356	1.356	1.572	1.264	0.976	1.219	0.984	0.992	0.877	0.867	0.867	1.529	3.316	3.21	3.21	1.912	2.08	2.08	1.274	1.084	1.028	2.231	2.336	2.336	0.922	0.312	0.869	0.869	0.644	0.865	1.432	0.23	0.226	0.26	0.26	1.078	1.038	0.958	0.958	1.644	1.925	2.007	2.432	2.915	3.088	3.325	2.522	2.434	1.964	1.237	1.146	1.843	0.768	0.653	0.467	0.946	1.239	1.154	1.154	0.832	0.488	0.175	0.214	0.21	0.11	0.106	0.172;
  10 675 1.952	0.924	4.404	4.448	4.196	4.228	3.748	4.332	2.108	2.324	2.996	0.98	0.988	3.9	3.4	3.408	2.132	0.292	0.388	0.264	0.268	0.316	0.264	0.284	0.284	0.236	0.324	0.172	0.228	0.348	0.168	0.228	0.252	0.252	0.244	0.496	0.308	0.3	0.296	0.208	0.3	0.208	0.236	0.224	0.184	0.148	0.148	0.216	0.34	0.336	0.372	0.34	0.564	0.564	1.32	0.808	0.981	0.968	0.76	0.548	0.631	0.468	0.388	0.412	0.532	0.532	3.4	1.568	0.488	0.488	0.372	0.528	1.844	3.668	3.436	3.192	2.024	2.1	2.052	1.572	3.156	0.948	4.268	3.592	2.168	0.896	0.992	0.992	0.988	0.356	0.356	0.28	0.32	0.336	0.436	0.436;
  11 680 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  12 684 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  13 692 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
]
#Reactive Demand
##Node# NodeName  T1	T2	T3	T4	T5	T6	T7	T8	T9	T10	T11	T12	T13	T14	T15	T16	T17	T18	T19	T20	T21	T22	T23	T24	T25	T26	T27	T28	T29	T30	T31	T32	T33	T34	T35	T36	T37	T38	T39	T40	T41	T42	T43	T44	T45	T46	T47	T48	T49	T50	T51	T52	T53	T54	T55	T56	T57	T58	T59	T60	T61	T62	T63	T64	T65	T66	T67	T68	T69	T70	T71	T72	T73	T74	T75	T76	T77	T78	T79	T80	T81	T82	T83	T84	T85	T86	T87	T88	T89	T90	T91	T92	T93	T94	T95	T96	T97	T98	T99	T100
Qit =[
  1  611 0.013014631	0	0	0.01859233	0.050199291	0	0	0.017352841	0	0.054537502	0	0	0	0	0	0	0.056396735	0.017972586	0	1.088890803	0	0	0.001859233	0.068171877	0.068171877	0	0	0	0	0	0.057636223	0.01859233	0	0	0	0	0	0.032226706	0.045861081	0	0	0	0	0	0	0.046480825	0.035325427	0	0.094201139	0.989731708	0	0	0	0.053298013	0.027888495	0	0	0	0	0	0.006197443	0.068791622	0.008056676	0	0	0	0	0.024789774	0.058255968	0	0	0	0	0	0	0	0.029747728	0	1.045508699	0	0	0	0.024789774	0.054537502	0	0	0	0	0	0.060115201	0.021691052	0	0	0	0	0;
  2  632 0.148738641	0.887473893	1.360958567	0.164851994	0.164851994	0.123329123	0.126427845	0.122089635	0.127667334	0.124568612	0.125188356	0.122089635	0.157415062	0.125188356	0.122709379	0.12146989	0.122089635	0.0656929	0.733777297	0.083665486	0.065073156	0.0656929	0.0656929	0.099159094	0.099159094	0.074369321	0.071270599	0.07003111	0.022310796	0.011775142	0.012394887	0.052058524	0.023550285	0.023550285	0.012394887	0.012394887	0.077468042	0.066932389	0.0656929	0.119610657	0.12146989	0.333422454	0.224967195	0.200797166	0.167330971	0.402214076	0.269588787	0.164851994	0.167950716	0.167330971	0.190881256	0.232404127	0.228685661	0.228065917	0.210713075	0.166091483	0.236742337	0.224347451	0.20575512	0.713945478	0.499513937	0.166091483	0.217530263	1.397523483	2.539092554	1.672070225	1.672070225	2.684732474	2.026563987	2.648787302	1.271715382	1.271715382	1.437806865	1.893938698	1.893938698	1.113060832	1.444004308	1.964589553	1.351042658	1.791061138	1.059762819	1.150245492	1.108102877	1.477470503	0.335901431	1.140329583	0.313590635	0.316069613	1.041790233	0.356972739	0.319788079	1.225854301	0.681718772	0.859585397	0.859585397	0.973618356;
  3  633 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  4  634 0.904826734	0.81000585	1.893318954	0.83851409	1.343605726	0.778398889	0.81992176	0.825499459	0.754848604	0.816823038	0.733777297	0.735016785	0.774680423	0.863923608	0.829217925	1.538205448	0.849049744	0.793272753	0.856486676	1.18990913	1.267377172	1.232051745	1.232051745	1.165739101	1.171936544	1.209740949	1.248784842	1.192388107	1.163879868	1.698099487	1.225854301	1.225854301	1.172556288	1.162020635	1.4285107	2.075523789	2.648787302	2.058170948	4.191950705	4.416298155	4.804258111	4.815413509	4.815413509	4.8711905	2.166626207	1.884642533	1.22709379	1.210980437	1.149006003	1.147766515	1.19858555	1.850556594	1.193627596	1.152724469	1.553079312	1.398143227	1.164499612	3.147681495	3.341661473	1.59274295	1.855514549	2.084819954	2.74422793	1.333689816	1.351042658	1.907573074	1.907573074	2.092256886	2.102172796	2.1455549	2.128202058	2.128202058	2.086059443	2.102172796	2.102172796	2.144935155	2.771496681	2.689690429	2.574417982	1.553699056	2.141836434	2.135019246	2.648787302	2.181500071	2.177781605	2.110849217	2.05135376	2.050114271	2.045156317	2.980970268	1.983801627	2.024085009	2.029042964	1.970786996	1.970786996	2.067467113;
  5  645 0.391058678	0.380523024	0.344577852	0.34705683	0.343338363	0.14997813	0.343338363	0.226206684	0.287561373	0.218150007	0.262151855	0.218150007	0.169809949	0.143780687	0.141301709	0.142541198	0.14006222	0.143780687	0.20575512	0.216290774	0.174148159	0.142541198	0.141301709	0.141301709	0.141301709	0.14006222	0.142541198	0.137583243	0.204515632	0.214431541	0.173528415	0.14006222	0.146259664	0.146259664	0.141921453	0.140681965	0.14006222	0.138822732	0.140681965	0.188402279	0.244179269	0.225586939	0.161133528	0.148738641	0.148738641	0.133245033	0.130146311	0.133864777	0.20575512	0.206994609	0.146259664	0.134484521	0.119610657	0.115272447	0.11279347	0.116511936	0.136963499	0.264630832	0.352014784	0.3228868	0.232404127	0.406552286	0.410270752	0.428243338	0.267729554	0.176627136	0.158654551	0.226826428	0.226826428	0.223107962	0.220628984	0.142541198	0.142541198	0.749270905	0.752369627	0.752369627	0.301195748	0.109075004	0.107835515	0.109075004	0.106596026	0.107835515	0.233643616	0.266490066	0.259053133	0.195839211	0.185303557	0.184064069	0.676760818	0.800709685	0.727579853	0.645773601	0.891192359	0.67428184	0.591236099	0.591236099;
  6  646 0.024170029	0.022930541	0.049579547	0.086144463	0.081186508	0.023550285	0.024789774	0.024170029	0.022310796	0.023550285	0.024789774	0.023550285	0.023550285	0.029747728	0.068791622	0.084904974	0.059495456	0.023550285	0.021691052	0.023550285	0.024789774	0.022310796	0.024789774	0.024170029	0.024170029	0.024789774	0.032226706	0.085524719	0.083665486	0.031606961	0.024789774	0.022310796	0.024789774	0.024789774	0.024789774	0.022310796	0.023550285	0.023550285	0.039663638	0.260292622	0.396636377	0.840373323	0.775919912	0.722002154	0.733157552	0.720762666	0.697212381	0.63709718	0.51562729	0.475963652	0.506331124	0.410890496	0.075608809	0.031606961	0.060734945	0.122089635	0.069411366	0.016733097	0.035325427	0.016733097	0.052678269	0.080566764	0.075608809	0.034705683	0.017352841	0.017352841	0.014873864	0.014873864	0.022930541	0.081186508	0.076848298	0.078087787	0.078087787	0.016113353	0.014873864	0.014873864	0.01859233	0.015493608	0.077468042	0.075608809	0.056396735	0.015493608	0.014873864	0.017352841	0.014873864	0.016113353	0.01425412	0.062594178	0.076848298	0.41708794	0.708367779	0.764144769	0.629660248	0.823020481	0.823020481	0.207614353;
  7  650 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  8  652 0.297477282	0.234883104	0.290660095	1.664013549	2.881191429	0.33714092	1.969547507	0.420806406	1.214698903	0.39849561	0.562108115	0.424524872	0.324746033	0.487738794	0.324126289	0.430102571	0.382382257	0.423905127	0.233643616	0.462949021	0.425144616	0.31421038	0.464188509	0.340239642	0.340239642	0.44993439	0.324746033	0.422665639	0.353874017	0.447455412	0.316689357	0.349535807	0.44993439	0.44993439	0.397875865	0.405312797	0.311731402	0.422665639	0.373086092	0.180965347	0.291899583	0.104117049	0.074369321	0.085524719	0.067552133	0.083045741	0.083045741	0.141301709	0.104117049	0.097299861	0.039663638	0.059495456	0.0656929	0.0656929	0.029127984	0.029127984	0.089243185	0.042762359	0.02850824	0.02850824	0.044001848	0.050199291	0.029127984	0.054537502	0.02850824	0.074369321	0.074369321	0.109694748	0.073749576	0.073749576	0.082425997	0.044621592	0.029747728	0.092961651	0.031606961	0.029747728	0.096060372	0.029127984	0.029747728	0.096680117	0.029127984	0.095440628	0.027888495	0.029747728	0.029127984	0.085524719	0.02850824	0.078087787	0.045241337	0.027888495	0.027888495	0.027888495	0.029127984	0.08862344	0.02850824	0.088003696;
  9  671 1.860472504	1.860472504	0.678000306	0.577601723	0.384861234	0.508190357	0.334042198	0.495175726	0.27144802	0.363789927	0.400354843	0.423285383	0.400354843	0.762285536	1.972026485	1.510936697	1.822048355	0.720762666	1.096327735	0.438778992	0.42576436	0.42576436	0.344577852	0.190881256	0.185303557	0.208853842	0.260292622	0.230544894	0.277645464	0.840373323	0.840373323	0.9742381	0.783356844	0.604870474	0.755468349	0.609828429	0.614786384	0.543515785	0.537318341	0.537318341	0.947589093	2.055072226	1.989379326	1.989379326	1.184951175	1.289068224	1.289068224	0.789554287	0.671802863	0.63709718	1.382649619	1.447722775	1.447722775	0.57140428	0.193360234	0.53855783	0.53855783	0.399115354	0.536078853	0.887473893	0.142541198	0.14006222	0.161133528	0.161133528	0.668084397	0.643294623	0.593715076	0.593715076	1.018859692	1.193007851	1.243826887	1.507218231	1.806554746	1.913770517	2.060649925	1.562995221	1.50845772	1.217177881	0.766623747	0.710227012	1.142188816	0.475963652	0.404693053	0.289420606	0.586278144	0.767863235	0.715184967	0.715184967	0.51562729	0.302435237	0.108455259	0.132625288	0.130146311	0.068171877	0.0656929	0.106596026;
  10 675 1.209740949	0.572643769	2.729354066	2.756622817	2.600447244	2.620279063	2.32280178	2.684732474	1.306421065	1.440285842	1.856754038	0.607349452	0.612307406	2.41700292	2.107130751	2.112088705	1.321294929	0.180965347	0.240460803	0.163612505	0.166091483	0.195839211	0.163612505	0.176007392	0.176007392	0.146259664	0.200797166	0.106596026	0.141301709	0.21567103	0.104117049	0.141301709	0.156175573	0.156175573	0.151217619	0.307393192	0.190881256	0.185923302	0.183444324	0.128906822	0.185923302	0.128906822	0.146259664	0.138822732	0.114032958	0.091722162	0.091722162	0.133864777	0.210713075	0.208234098	0.230544894	0.210713075	0.349535807	0.349535807	0.818062527	0.500753425	0.607969196	0.59991252	0.471005697	0.339619897	0.391058678	0.29004035	0.240460803	0.255334667	0.329703988	0.329703988	2.107130751	0.971759123	0.302435237	0.302435237	0.230544894	0.327225011	1.14280856	2.273222233	2.129441547	1.978223928	1.254362541	1.301463111	1.271715382	0.9742381	1.955913132	0.587517633	2.645068836	2.226121664	1.343605726	0.555290927	0.614786384	0.614786384	0.612307406	0.220628984	0.220628984	0.173528415	0.198318188	0.208234098	0.270208532	0.270208532;
  11 680 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  12 684 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
  13 692 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
]

#%Prated1 include values >=0 and  <= Prated1
#%Prated2 include values >Prated1 and  <= Prated2
#CT	  Type_indicator   CTDO   CTCO	  CTB   FC     VC     OC2    OC1    OC0    Pmax    Pmin   Δc   StandbyLoss(L0) %Prated1 eff1 %Prated2 eff2
CTSet =[
     1	    1            0	      0	    1	    10000.0  300.0	10.0   5.0    2.0    100.0    -100.0   10.0     -0.3         0.50    0.95   1.0     0.90;
     2	    1            0	      1	    0	    20000.0	 250.0	20.0   10.0   4.0    100.0    0.0   10.0     -0.0         0.50    0.95   1.0     0.90;
     3	    2            0	      1	    0	    25000.0	 200.0	30.0   15.0   8.0    100.0    0.0   10.0     -0.0         0.50    0.95   1.0     0.90;
     4	    1            0	      1	    0	    30000.0	 150.0	40.0   20.0   10.0   100.0    0.0   10.0     -0.0         0.50    0.95   1.0     0.90;
     5	    2            0	      1	    0	    35000.0	 100.0	50.0   25.0   5.0    100.0    0.0   10.0     -0.0         0.50    0.95   1.0     0.90;
]

   #DT   DTDO   DTCO   DTB   FC     VC     OC2    OC1    OC0     UT   DT   RampUpR   RampDownR   MaxP   MaxQ   MinP   MinQ   Δd StandbyLoss(L0) %Prated1 eff1 %Prated2 eff2
DTSet=[
    1	    1   	0	     0	  20000.0	  0.0	  50.0    25.0    6.0    3	  3      200	    200	      250	    250	  -0	-250    10.0  -0.0         0.50    0.95   1.0     0.90;
    2	    0   	1	     0	  10000.0	  0.0	  40.0    20.0    5.0    0    0      0        0         275	    250	  -0	-250    10.0  -0.0         0.50    0.95   1.0     0.90;
    3	    0   	1	     0	  25000.0	  0.0	  30.0    15.0    4.0    0    0      0        0         300	    250	  -0	-250    10.0  -0.0         0.50    0.95   1.0     0.90;
    4	    1   	0	     0	  30000.0	  0.0	  20.0    10.0    3.0    3	  3      200	    200	      225	    250	  -0	-250    10.0  -0.0         0.50    0.95   1.0     0.90;
    5	    1   	0	     0	  35000.0	  0.0	  10.0    5.0    	2.0    3    3      200	    200       200	    250	  -0	-250    10.0  -0.0         0.50    0.95   1.0     0.90;
]

  #ArcID   Start   End   StartNodeId   EndnodeID   Length     Resistance     Reactance                    Tij 	newline cost
Linedet =[
    1	   650   632   7    2	   0.37878800   0.129469738   0.391477398	  200	0	0   ;
    2	   632   633   2    3    0.09469700   0.070823886   0.113349152	  200	0	0   ;
    3	   632   645   2    5	   0.09469700   0.125625040   0.128030344	  200	0	0   ;
    4	   632   671   2    9    0.37878800   0.129469738   0.391477398	  200	0	0   ;
    5	   645   646   5    6    0.05681820   0.075375024   0.076818206	  200	0	0   ;
    6	   633   634   3    4    0.00189394   0.000647349   0.001957387	  200	0	0   ;
    7	   671   692   9    13   0.00189394   0.000647349   0.001957387	  200	0	0   ;
    8	   671   680   9    11   0.18939400   0.064734869   0.195738699	  200	0	0   ;
    9	   671   684   9    12   0.05681820   0.075375024   0.076818206	  200	0	0   ;
    10   692   675   13   10   0.09469700   0.075299898   0.040931200   200	0	0   ;
    11   684   652   12   8    0.15151520   0.203409156   0.077636388	  200	0	0   ;
    12   684   611   12   1    0.05681820   0.037761376   0.038281262	  200	0	0   ;
    13   646   611   6    1    0.37878800   0.129469738   0.391477398   200 1 1000 ;
    14	 650   632   7    2	   0.37878800   0.129469738   0.391477398	  200	1 1000 ;
    15	 632   633   2    3    0.09469700   0.070823886   0.113349152	  200	1 1000 ;
    16	 632   645   2    5	   0.09469700   0.125625040   0.128030344	  200	1 1000 ;
    17	 632   671   2    9    0.37878800   0.129469738   0.391477398	  200	1 1000 ;
    18	 645   646   5    6    0.05681820   0.075375024   0.076818206	  200	1 1000 ;
    19	 633   634   3    4    0.00189394   0.000647349   0.001957387	  200	1 1000;
    20	 671   692   9    13   0.00189394   0.000647349   0.001957387	  200	1 1000 ;
    21	 671   680   9    11   0.18939400   0.064734869   0.195738699	  200	1 1000 ;
    22	 671   684   9    12   0.05681820   0.075375024   0.076818206	  200	1 1000 ;
    23   692   675   13   10   0.09469700   0.075299898   0.040931200   200	1 1000 ;
    24   684   652   12   8    0.15151520   0.203409156   0.077636388	  200	1 1000 ;
    25   684   611   12   1    0.05681820   0.037761376   0.038281262	  200	1 1000 ;
]
	#Existing Continuous generation
  #Node ContinuousID-->
Bgc_ex = []

	#Existing Discrete generation
  #Node Slot#    DiscreteID-->   D1   D2   D3  D4  D5
Bgd_ex = [
#  1	    1 0 0 0 0 0
#  2 	    1 0 0 0 0 0
#  3 	    1 0 0 0 0 0
#  4      1 0 0 0 0 0
#  5      1 0 1 0 0 0
#  5      2 0 0 0 0 0
#  5      3 0 0 0 0 0
#  6      1 0 0 0 0 0
#  7      1 0 0 0 0 0
#  8      1 0 1 0 0 0
#  9      1 0 0 0 0 0
#  10     1 0 0 0 0 0
#  11     1 0 0 0 0 0
#  12     1 0 0 0 0 0
#  13     1 0 0 0 0 0
  ]
