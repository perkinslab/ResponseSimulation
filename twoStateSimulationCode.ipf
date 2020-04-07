#pragma rtGlobals=3		

///////////
//Annotated simulation code to accompany:
//
//"Correcting molecular transition rates measured by single-molecule force 
//spectroscopy for limited temporal resolution"
//
//by
//
//David R. Jacobson and Thomas T. Perkins (2020)
//////////

//////////
//INSTRUCTIONS:
//This code simulates an atomic force microscopy experiment in which a biomolecule is repeatedly 
//stretched and relaxed.  The molecule may exist in either of two states (i.e., folded and unfolded). 
//Transition rates between the states vary as a function of force and are given by the Dudko rate equation
//(Dudko et al., Phys. Rev. Lett. 96, 108101 [2006]), which is parameterized in terms of the zero-force 
//transition rate, transition barrier height, and distance to the transition state of the one-dimensional free-
//energy landscape separating the states.  After simulating force ramps until a specified number of state
//transitions are recorded, the effect of a non-zero instrument time resolution is modeled by eliminating 
//transitions corresponding to dwells in states shorter than that time.
//
//Running the simulation code consists of  five main steps, each a separate function call:
//	(1) setParams(), which creates a table in which the user specifies the simulation parameters (e.g., 
//		free-energy landscape parameters, AFM pulling parameters, number of iterations, response times
//		to model)
//	(2) simBatch(), which runs the simulation.  This step is slow.  For example, on a typical PC, a simulation 
//		of 100,000 transitions in each direction and 7 time response values will take ~27 hours.
//	(3) calcRateMap(bootstraps), which calculates the force-dependent transition rates from the simulated
//		data.  Bootstrapping is used to estimate errors, with the bootstraps parameter  being the number of
//		bootstrap iterations to use (typically ~50).
//	(4) correctRateMap(), which applies the time-response corrections developed in our paper.
//	(5) makePlots(filter), which produces plots like Figs. 2(b), 2(c), 5(a), 5(b), 5(c), and 6 of the paper. The 
//		parameter filter is an index of the filters wave that specifies which time-response value to use in 
//		Figs. 5 and 6.
//
//Detailed comments describing how to use each function are given alongside each function below.  These
//five key functions are given below in the order in which they are to be called, followed but other ancillary 
//functions
//////////////


Function setParams()
//This function makes two tables of simulation parameters and sets the parameters to the default values
//we used in the paper.  If the simulations are to be run with different parameters, those should be entered
//into these tables before calling the next function.

//The filters wave is a list of instrument response times to be modeled, each specified in microseconds.  The 
//response time specified as 0 in this list will correspond to modeling a time response equal to the simulation 
//time step (0.1 us by default).

//The Param wave contains the other simulation parameters, each designated by the corresponding label
//in the ParamLabels text wave.  These parameters are:
//
//	0	Timestep (s)		Discrete simulation time step, which should be fast compared to transition rates
//	1	Num. transitions	Simulations will run until there are this many transitions recorded in each direction
//	2	Bin size (N)		Force bin size for calculating rates
//	3	Velocity (m/s)	Rate at which simulated AFM cantilever is moved away from/toward surface
//	4	Start Z (m)		Closest approach of AFM cantilever to surface
//	5	End Z (m)		Furthest distance of AFM cantilever from surface
//	6	Spring const. (N/m) Spring constant of AFM cantilever
//	7	Linker Lc (m)		Contour length of the molecule in folded state
//	8	Liberated Lc (m) 	Contour length *increase* upon unfolding
//	9	dx_01 (m)		Distance to transition state - unfolding
//	10	dx_10 (m) 		Distance to transition state - refolding
//	11	k0_01 (1/s)		Zero-force transition rate - unfolding
//	12	k0_10 (1/s)		Zero-force transition rate - refolding
//	13	dG_01 (J)		Transition barrier height - unfolding
//	14 	dG_10 (J)		Transition barrier height - refolding

//Several ancillary functions can be used to help choose parameters:
//
//	calcDudko(plot=1) makes a plot of the expected force dependence of the rates for the current parameters
//
//	simSingle(1,plot=1) shows a single simulated unfolding force-extension curve with the current parameters
	
	//Create parameter wave
	Make/O/T ParamLabels = {"Timestep (s)", "Num. transitions", "Bin size (N)", "Velocity (m/s)", "Start Z (m)", "End Z (m)", "Spring const. (N/m)", "Linker Lc (m)", "Liberated Lc (m)", "dx_01 (m)", "dx_10 (m)", "k0_01 (1/s)", "k0_10 (1/s)", "dG_01 (J)", "dG_10 (J)"}
	Make/O Param = {1e-7, 1000, 10e-12, 300e-9, 5e-9, 2e-8, 20e-3, 12.6e-9, 1e-9, 0.4e-9, -0.2e-9, 20, 300000, 4.4e-20, 4e-21}
	
	//Create wave of filter times
	Make/O/T filterLabel = {"Filtering times (us)"}
	Make/O filters = {0,10,20,50,100,200,500,1000}
	
	//Display both tables
	Edit ParamLabels, Param
	Edit filterLabel, filters
End

Function simBatch()
//This function continues to simulate force ramps of increasing and decreasing force until the number of transitions 
//in each direction reaches the user-specified amount.  The following quantities are 
//recorded (1) The amount of time spent in each state in each force bin, (2) the number of transitions in each 
//force bin, (3) the number of trajectories in the state at a particular force at the center of the force bin; i.e.,
//N_i in Eq. 3 of the paper.

	//Reference user-specified parameter waves:
	Wave param, filters
	
	Variable transitions_desired = param[1]
	
	//Set certain simulation parameters
	Variable curr_direction = 1
	Variable curr_ramp = 0
	Variable num_trans_sofar = 0
	Variable num_force_bins = ceil( (Param[5]-calcXfromZ(Param[5], Param[7]))*Param[6] / Param[2]) //Approximately 
	Variable binsize = Param[2]
	
	//Make the waves that will store the values calculated, separately, for each filtering time. Each wave will be
	//three-dimensional, with these dimensions: 0- force bin, 1- record, 2- filter
	Make/O/N=(num_force_bins,0,numpnts(filters)) times_0, times_1, trans_01, trans_10, N_0,  N_1
	Make/O/N=0 ramp_direction

	//Iteration over single force-ramp simulations, repated until the desired number of transitions are reached
	Variable k, n, currbin
	do
		//Run a single simulation
		simSingle(curr_direction)
		Wave F_sim,  state_sim //simulation output
		
		//Add a new entry in the appropriate dimension ot the 3D waves:
		Redimension/N=(num_force_bins, curr_ramp+1, numpnts(filters))   times_0, times_1, trans_01, trans_10, N_0,  N_1, ramp_direction
		Redimension/N=(curr_ramp+1) ramp_direction
		ramp_direction[curr_ramp] = curr_direction
		
		//Parse the results separately for each filtering time
		for (k=0; k<numpnts(filters); k+=1) //filters indexed by k
			if (filters[k] == 0) //i.e., unfiltered
				Duplicate/O F_sim, F_sim_ds
				Duplicate/O state_sim, state_sim_ds
			else
				filterRecord(F_sim, state_sim, filters[k])
			endif
			//Make the temporary waves to store the results of this one filtering.  Each is a 1D wave,
			//with each entry a different force value
			Make/O/N=(num_force_bins) curr_times_0, curr_times_1, curr_trans_01, curr_trans_10, curr_N_0,  curr_N_1
			curr_times_0 = 0
			curr_times_1 = 0
			curr_trans_01 = 0
			curr_trans_10 = 0
			curr_N_0 = 0
			curr_N_1 = 0
			//Iterate through the simulation wave and find the dwells and transitions
			for (n=0; n<numpnts(F_sim_ds)-1; n+=1)
				currbin = floor(F_sim_ds[n]/binsize)
				//Assign dwells
				if (state_sim_ds[n] == 1)
					curr_times_1[currbin] += 1
				else
					curr_times_0[currbin] += 1 
				endif
				//Assign transitions
				if (state_sim_ds[n] == 1 && state_sim_ds[n+1] == 0)
					curr_trans_10[currbin] += 1
				elseif (state_sim_ds[n] == 0 && state_sim_ds[n+1] == 1)
					curr_trans_01[currbin] += 1
				endif
				//Figure if this point contributes to N_i (i.e., if, between this point and the next, the 
				//dwell crosses the critical force associated with the bin)
				if (state_sim_ds[n] == state_sim_ds[n+1]) //i.e., not a transition
					//Is this a crossing of a halfway force?
					if ( abs( floor( (F_sim_ds[n] - (binsize/2))/binsize) - floor( (F_sim_ds[n+1] - (binsize/2))/binsize)) == 1)
						if (state_sim_ds[n] == 1)
							curr_N_1[currbin] += 1
						else
							curr_N_0[currbin] += 1
						endif
					endif
				endif
			endfor
			
			//Increment the total number of transitions only with the lowest-filter result (assumed to be first)
			if (k == 0)
				num_trans_sofar += sum(curr_trans_01)
			endif 
			
			//Save the resulting values to the correct places in the 3D waves
			for (n=0; n<num_force_bins; n+=1)
				times_0[n][curr_ramp][k] = curr_times_0[n]
				times_1[n][curr_ramp][k] = curr_times_1[n]
				trans_01[n][curr_ramp][k] = curr_trans_01[n]
				trans_10[n][curr_ramp][k] = curr_trans_10[n]
				N_0[n][curr_ramp][k] = curr_N_0[n]
				N_1[n][curr_ramp][k] = curr_N_1[n]
			endfor
			
		endfor
		
		//Set values for next cycle
		curr_direction = curr_direction*(-1)
		curr_ramp += 1
		
		//Optional: Save as you go: (Because simBatch() takes a long time to run and IGOR could crash)
		SaveExperiment
		Print num_trans_sofar
		
	while (num_trans_sofar < transitions_desired)
	
	//Convert the times from timestep increments into seconds:
	times_0 = times_0*Param[0]
	times_1 = times_1*Param[0]
	
	//Cleanup:
	KillWaves curr_times_0, curr_times_1, curr_trans_01, curr_trans_10, curr_N_0,  curr_N_1
End


Function calcRateMap(bootstraps)
//After running simBatch to generate the 3D waves containing information about the transitions, this function 
//processes them to get the rate maps.  The rates are calculated both using the Zhang-Dudko method (Eq. 3 of
//paper) and as the number of transitions per time (Eq. 4 of paper).  Bootstrapping is done on a per-rate-map 
//basis to get error bars from both methods. The bootstraps parameter specifies the number of samplings to do in
//the bootstrap analysis; 50 is a suggested value.
	Variable bootstraps
	
	Wave times_0, times_1, trans_01, trans_10, N_0, N_1 //outputs of simulation
													 // Dimensions:
													 //0- force bin, 1- ramp, 2- filter
	//User-specified simulation parameters
	Wave Param, filters
	
	//Calculate certain additional parameters
	Variable num_ramps = dimSize(times_0, 1)
	Variable num_force_bins = dimSize(times_0, 0)
	Variable num_filters = dimSize(times_0,2)
	
	//Make waves that will store the rate maps:
	Make/O/N=(num_force_bins, num_filters) k_01_ZD, k_10_ZD, k_01_J, k_10_J, dk_01_ZD, dk_10_ZD, dk_01_J, dk_10_J
	SetScale/P x, (Param[2]/2), Param[2], k_01_ZD, k_10_ZD, k_01_J, k_10_J, dk_01_ZD, dk_10_ZD, dk_01_J, dk_10_J
	
	//Iterate through the data to construct the maps:
	Make/O/N=(num_ramps) bt_index
	Variable ft, bt, k, n
	for (ft=0; ft<num_filters; ft+=1) //iteration over different time-response filters
		//Make waves to store the rates from each iteration of the bootstrapping
		Make/O/N=(num_force_bins, bootstraps) k_01_ZD_curr, k_10_ZD_curr, k_01_J_curr, k_10_J_curr
		for (bt=0; bt<bootstraps; bt+=1) //iteration over bootstraps
			//Populate the random bootstrap indecies 
			for (k=0; k<numpnts(bt_index); k+=1)
				bt_index[k] = floor(abs(enoise(num_ramps)))
			endfor
			
			//Add up cumulative values across these chosen indecies
			Make/O/N=(num_force_bins) times_0_cum, times_1_cum, trans_01_cum, trans_10_cum, N_0_cum, N_1_cum
			times_0_cum = 0
			times_1_cum = 0
			trans_01_cum = 0
			trans_10_cum = 0
			N_0_cum = 0
			N_1_cum = 0
			for (k=0; k<numpnts(bt_index); k+=1) //over ramps
				for (n=0; n<num_force_bins; n+=1) //over force bins
					times_0_cum[n] += times_0[n][bt_index[k]][ft]
					times_1_cum[n] += times_1[n][bt_index[k]][ft]
					trans_01_cum[n] += trans_01[n][bt_index[k]][ft]
					trans_10_cum[n] += trans_10[n][bt_index[k]][ft]
					N_0_cum[n] += N_0[n][bt_index[k]][ft]
					N_1_cum[n] += N_1[n][bt_index[k]][ft]
				endfor
			endfor
			
			//Calculate the rate both ways from this particular bootstrapping
			for (k=0; k<num_force_bins; k+=1)
				//Metod of this paper (Eq. 4)
				k_01_J_curr[k][bt] = trans_01_cum[k]/times_0_cum[k]
				k_10_J_curr[k][bt] = trans_10_cum[k]/times_1_cum[k]
				
				//Zhang-Dudko method (Eq. 3)
				k_01_ZD_curr[k][bt] = estLoadingRate(Param[2]*(1/2+k), Param[7])*trans_01_cum[k]/N_0_cum[k]/Param[2]
				k_10_ZD_curr[k][bt] = estLoadingRate(Param[2]*(1/2+k), Param[7]+Param[8])*trans_10_cum[k]/N_1_cum[k]/Param[2]

			endfor
			
			//Cleanup:
			KillWaves times_0_cum, times_1_cum, trans_01_cum, trans_10_cum, N_0_cum, N_1_cum
		endfor
		
		//Calculate the values and standard errors:
		Make/O/N=(bootstraps) k_01_J_temp, k_10_J_temp, k_01_ZD_temp, k_10_ZD_temp
		for (k=0; k<num_force_bins; k+=1)
			k_01_J_temp = 0
			k_10_J_temp =0
			k_01_ZD_temp = 0
			k_10_ZD_temp = 0
			for (n=0; n<bootstraps; n+=1)
				k_01_J_temp[n] = k_01_J_curr[k][n]
				k_10_J_temp[n] = k_10_J_curr[k][n]
				k_01_ZD_temp[n] = k_01_ZD_curr[k][n]
				k_10_ZD_temp[n] = k_10_ZD_curr[k][n]
			endfor
			k_01_J[k][ft] = mean(k_01_J_temp)
			dk_01_J[k][ft] = sqrt(variance(k_01_J_temp))
			k_10_J[k][ft] = mean(k_10_J_temp)
			dk_10_J[k][ft] = sqrt(variance(k_10_J_temp))
			
			k_01_ZD[k][ft] = mean(k_01_ZD_temp)
			dk_01_ZD[k][ft] = sqrt(variance(k_01_ZD_temp))
			k_10_ZD[k][ft] = mean(k_10_ZD_temp)
			dk_10_ZD[k][ft] = sqrt(variance(k_10_ZD_temp))
		endfor
		
		//Cleanup:
		KillWaves k_01_J_temp, k_10_J_temp, k_01_ZD_temp, k_10_ZD_temp, k_01_ZD_curr, k_10_ZD_curr, k_01_J_curr, k_10_J_curr
	endfor
End

Function correctRateMap()
//Takes the output of calcRateMap() and applies the two corrections developed in the paper.  Generates 
//waves  giving (1) the effect of the high-k correction (Eq. 7) and (2) the effects of both corrections (Eqs.
//7 and 8).
//Note: This code assumes that the first filter is the zero filter
	
	//Load simulation parameters, rate maps, and dwell times in states.  The correction will be applied to
	//rates calculated using Eq. 4.
	Wave Param, filters
	Wave k_01_J, k_10_J 
	Wave times_0, times_1 //(from simBatch())
	
	//Calculate the number of transitions occuring in each rate bin from the rate and the time spent in the state:
	Variable ft, k, bn, n
	Variable T0, T1
	Duplicate/O k_01_J, N_01_J
	Duplicate/O k_10_J, N_10_J
	for (ft=0; ft<numpnts(filters); ft+=1)
		for (bn=0; bn<dimSize(k_01_J, 0); bn+=1)
			T0 = 0
			T1 = 0
			for (n=0; n<dimSize(times_0, 1); n+=1)
				T0 += times_0[bn][n][ft]
				T1 += times_1[bn][n][ft]
			endfor
			N_01_J[bn][ft] = k_01_J[bn][ft]*T0
			N_10_J[bn][ft] = k_10_J[bn][ft]*T1
		endfor
	endfor
	
	//HIGH-K CORRECTION (Eq. 7)
	
	//Create waves to store rates after this first correction is applied
	Duplicate/O k_01_J, k_01_J_HR
	Duplicate/O k_10_J, k_10_J_HR //"HR" is the "high rate" correction
	
	//Apply the correction for all filter times and rates:
	for (ft=1; ft<numpnts(filters); ft+=1)
		for (k=0; k<dimSize(k_01_J_HR,0); k+=1)
			k_01_J_HR[k][ft] = k_01_J[k][ft]/(1-k_01_J[k][ft]*filters[ft]*1e-6)
			k_10_J_HR[k][ft] = k_10_J[k][ft]/(1-k_10_J[k][ft]*filters[ft]*1e-6)
		endfor
	endfor
	
	//CORRESPONDING LOW-K CORRECTION (Eq. 8)
	
	//Find the number of missing transitions that were added to each bin by the first correction
	
	//Create additional waves to store results and intermediate quantities:
	Duplicate/O k_01_J_HR, N_missing_01_J, N_toadd_01_J, k_01_J_fullcorr
	Duplicate/O k_10_J_HR, N_missing_10_J, N_toadd_10_J, k_10_J_fullcorr
	
	//Initialize:
	N_toadd_10_J = 0
	N_toadd_01_J = 0
	
	//Iterate through filter times and rates.  For each, find the number of missing transitions implied 
	//by the first correction:
	for (ft=1; ft<numpnts(filters); ft+=1)
		for (k=0; k<dimSize(N_missing_01_J, 0); k+=1)
			T0 = 0
			T1 = 0
			for (n=0; n<dimSize(times_0, 1); n+=1)
				T0 += times_0[k][n][ft]
				T1 += times_1[k][n][ft]
			endfor
			N_missing_01_J[k][ft] = k_01_J_HR[k][ft]/k_01_J[k][ft]*(k_01_J_HR[k][ft]-k_01_J[k][ft])*T0
			N_missing_10_J[k][ft] = k_10_J_HR[k][ft]/k_10_J[k][ft]*(k_10_J_HR[k][ft]-k_10_J[k][ft])*T1
		endfor
	endfor
	
	//Assign these missing N to the correct bins in the other trace, accounting for the force drop at unfolding
	//(which varies as a function of force)
	
	//Evaluate the wormlike chain elasticity of the folded (f) and unfolded states (u)
	Make/O/N=(10000) WLC_u, WLC_f
	SetScale/I x, Param[4], Param[5], WLC_u, WLC_f
	WLC_u = calcXfromZ(x, Param[7]+Param[8])
	WLC_u = (x-WLC_u)*Param[6]
	WLC_f = calcXfromZ(x, Param[7])
	WLC_f = (x-WLC_f)*Param[6]
	
	//Iterate over force bins and, for each one, find which bins in the WLC of the other state it maps to:
	Variable bin_start_F, bin_end_F, reg_start_pnt, reg_end_pnt, reg_length, curr_bin, addvalue
	Make/O/N=(dimSize(N_missing_10_J, 0)) distribution 
	for (bn= 0; bn<dimSize(N_missing_01_J, 0); bn+=1)
		bin_start_F = bn*Param[2]
		bin_end_F = (bn+1)*Param[2]
		distribution = 0
		//Find the point numbers in the WLC corresponding to the start and end of the bin region
		reg_start_pnt = 0
		reg_end_pnt = numpnts(WLC_u) - 1 //default values
		for (k=0; k<numpnts(WLC_u); k+=1)
			if (WLC_u[k] >= bin_start_F)
				reg_start_pnt = k
				break
			endif
		endfor
		for (k=reg_start_pnt; k<numpnts(WLC_u); k+=1)
			if (WLC_u[k] >= bin_end_F)
				reg_end_pnt = k
				break
			endif
		endfor
		reg_length = reg_end_pnt - reg_start_pnt
		for (k=reg_start_pnt; k<= reg_end_pnt; k+=1)
			curr_bin = floor(WLC_f[k]/Param[2])
			distribution[curr_bin] += 1
		endfor
		distribution = distribution/reg_length
		for (ft=1; ft<numpnts(filters); ft+=1)
			for (k=0; k<dimsize(N_toadd_01_J, 0); k+=1)
				addvalue = (distribution[k]*N_missing_10_J[bn][ft])
				if (numtype(addvalue) == 2)
					addvalue = 0
				endif
				N_toadd_01_J[k][ft] = N_toadd_01_J[k][ft] + addvalue
			endfor
		endfor
	endfor
	//Repeat for the other transition
	for (bn= 0; bn<dimSize(N_missing_01_J, 0); bn+=1)
		bin_start_F = bn*Param[2]
		bin_end_F = (bn+1)*Param[2]
		distribution = 0
		//Find the point numbers in the WLC corresponding to the start and end of the bin region
		reg_start_pnt = 0
		reg_end_pnt = numpnts(WLC_f) - 1 //default values
		for (k=0; k<numpnts(WLC_f); k+=1)
			if (WLC_f[k] >= bin_start_F)
				reg_start_pnt = k
				break
			endif
		endfor
		for (k=reg_start_pnt; k<numpnts(WLC_f); k+=1)
			if (WLC_f[k] >= bin_end_F)
				reg_end_pnt = k
				break
			endif
		endfor
		reg_length = reg_end_pnt - reg_start_pnt
		for (k=reg_start_pnt; k<= reg_end_pnt; k+=1)
			curr_bin = floor(WLC_u[k]/Param[2])
			distribution[curr_bin] += 1
		endfor
		distribution = distribution/reg_length
		for (ft=1; ft<numpnts(filters); ft+=1)
			for (k=0; k<dimsize(N_toadd_10_J, 0); k+=1)
				addvalue = (distribution[k]*N_missing_01_J[bn][ft])
				if (numtype(addvalue) == 2)
					addvalue = 0
				endif
				N_toadd_10_J[k][ft] = N_toadd_10_J[k][ft] + addvalue
			endfor
		endfor
	endfor
	
	//Add these toadd values (the number of missed transition to be added to a certain bin in this correction)
	//to the relevant rates to get the low-k corrected rates:
	for (ft=1; ft<numpnts(filters); ft+=1)
		for (bn=0; bn<dimSize(k_01_J_fullcorr, 0); bn+=1)
			T0 = 0
			T1 = 0
			for (n=0; n<dimSize(times_0, 1); n+=1)
				T0 += times_0[bn][n][ft]
				T1 += times_1[bn][n][ft]
			endfor
			k_01_J_fullcorr[bn][ft] = k_01_J_fullcorr[bn][ft] + N_toadd_01_J[bn][ft]/T0
			k_10_J_fullcorr[bn][ft] = k_10_J_fullcorr[bn][ft] + N_toadd_10_J[bn][ft]/T1
		endfor
	endfor
	
	//Cleanup:
	KillWaves WLC_u, WLC_f, distribution 
End

Function makePlots(filter)
//Produces plots like Figs. 2(b), 2(c), 5(a), 5(b), 5(c), and 6 of the paper. The filter parameter is an 
//index of the filters wave that specifies which time-response value to use in Figs. 5 and 6.
	Variable filter
	
	//Simulation and correction outputs:
	Wave filters
	Wave k_01_J, k_10_J
	Wave dk_01_J, dk_10_J
	Wave k_01_ZD, k_10_ZD
	Wave dk_01_ZD, dk_10_ZD4
	Wave k_01_J_HR, k_10_J_HR
	Wave k_01_J_fullcorr, k_10_J_fullcorr
	
	//Figure 2(b): Comparison of rates from simulation results without filtering, calculated using Eq. 3
	//or Eq. 4, with the theoretically expected rate maps:
		
		//Calculate the theoretical curves
		calcDudko()
		Wave k_01_ideal, k_10_ideal //output
		
		//Add these to the plot
		Display/N=Fig2b k_01_ideal
		ModifyGraph log(left)=1
		Appendtograph k_10_ideal
		ModifyGraph rgb=(34816,34816,34816)
		
		//Add calculated rates to the plot
		Appendtograph k_01_J[][0]
		Appendtograph k_10_J[][0]
		Appendtograph k_01_ZD[][0]
		Appendtograph k_10_ZD[][0]
		ModifyGraph mode(k_01_J)=3,marker(k_01_J)=8,mode(k_10_J)=3,marker(k_10_J)=8;DelayUpdate
		ModifyGraph rgb(k_10_J)=(0,0,65280),mode(k_01_ZD)=3,marker(k_01_ZD)=1;DelayUpdate
		ModifyGraph mode(k_10_ZD)=3,marker(k_10_ZD)=1,rgb(k_10_ZD)=(0,0,65280)
		
		//Add labels
		Label left "k (1/s)";DelayUpdate
		Label bottom "Force (N)"
		Legend/C/N=text0/J/F=0/A=MC "\\s(k_01_ideal) Expected rates\r\\s(k_01_ZD) Unfolding (Eq. 3)\r\\s(k_10_ZD) Refolding (Eq. 3)\r\\s(k_01_J) Unfolding (Eq. 4)";DelayUpdate
		AppendText "\\s(k_10_J) Refolding (Eq. 4)"
		
	//Figure 2(c): Plot of simulated rate maps (calculated using Eq. 4) for all response-time values:
	
		//Iterate through all filter values and plot each:
		Variable ft
		Display/N=Fig2c
		for (ft=0; ft<numpnts(filters); ft+=1)
			//Plot rates
			Appendtograph k_01_J[][ft]/TN=$("unfold"+num2str(ft))
			Appendtograph k_10_J[][ft]/TN=$("refold"+num2str(ft))
			ModifyGraph rgb($("refold"+num2str(ft)))=(0,0,65280)
			//Add error bars
			Errorbars $("unfold"+num2str(ft)), Y wave=(dk_01_J[][ft], dk_01_J[][ft])
			Errorbars $("refold"+num2str(ft)), Y wave=(dk_10_J[][ft], dk_10_J[][ft])
		endfor
		
		//Add labels
		ModifyGraph log(left)=1
		Label left "k_obs (1/s)";DelayUpdate
		Label bottom "Force (N)"
		
	//Figure 5(a): Comparison of time-response-biased rate map to expected values.  Shows data for the 
	//particular time-response value corresponding to the entry in the wave filters with index filter
	
		//Plot expected curves
		Display/N=Fig5a k_01_ideal
		ModifyGraph log(left)=1
		Appendtograph k_10_ideal
		ModifyGraph rgb=(34816,34816,34816)
		
		//Append biased rates
		Appendtograph k_01_J[][filter]
		Appendtograph k_10_J[][filter]
		ModifyGraph mode(k_01_J)=3,marker(k_01_J)=19,mode(k_10_J)=3,marker(k_10_J)=19;DelayUpdate
		ModifyGraph rgb(k_10_J)=(0,0,65280)
		
		//Add labels
		Label left "k (1/s)";DelayUpdate
		Label bottom "Force (N)"
		
	//Figure 5(b): Same as Fig. 5(a), but now plotting the data after the first correction (Eq. 7) is applied
	
		//Plot expected curves
		Display/N=Fig5b k_01_ideal
		ModifyGraph log(left)=1
		Appendtograph k_10_ideal
		ModifyGraph rgb=(34816,34816,34816)
		
		//Append biased rates
		Appendtograph k_01_J_HR[][filter]
		Appendtograph k_10_J_HR[][filter]
		ModifyGraph mode(k_01_J_HR)=3,marker(k_01_J_HR)=16,mode(k_10_J_HR)=3;DelayUpdate
		ModifyGraph marker(k_10_J_HR)=16,rgb(k_10_J_HR)=(0,0,65280)
		
		//Add labels
		Label left "k (1/s)";DelayUpdate
		Label bottom "Force (N)"
		
	//Figure 5(c): Same as Fig. 5(a), but now plotting the data after both corrections are applied
	
		//Plot expected curves
		Display/N=Fig5c k_01_ideal
		ModifyGraph log(left)=1
		Appendtograph k_10_ideal
		ModifyGraph rgb=(34816,34816,34816)
		
		//Append biased rates
		Appendtograph k_01_J_fullcorr[][filter]
		Appendtograph k_10_J_fullcorr[][filter]
		ModifyGraph mode(k_01_J_fullcorr)=3,marker(k_01_J_fullcorr)=18,mode(k_10_J_fullcorr)=3;DelayUpdate
		ModifyGraph marker(k_10_J_fullcorr)=18,rgb(k_10_J_fullcorr)=(0,0,65280)
		
		//Add labels
		Label left "k (1/s)";DelayUpdate
		Label bottom "Force (N)"
		
	//Figure 6: Plot of the correction factor for each rate for the response time of interest.  The correction
	//factor is the ratio of the number of points added to each bin by the second correction (Eq. 8) to the 
	//number of points in that bin to begin with.
		
		//Load waves relevant to this calculation
		Wave N_01_J, N_10_J
		Wave N_toadd_01_J, N_toadd_10_J
	
		//Make waves to store result and then calculate the correction factor
		Make/N=(dimsize(N_01_J,0))/O CF_01_J
		Make/N=(dimsize(N_10_J,0))/O CF_10_J
		SetScale/P x, dimOffset(N_01_J,0), dimDelta(N_01_J, 0), CF_01_J, CF_10_J
		Variable k
		for (k=0; k<numpnts(CF_01_J); k+=1)
			CF_01_J[k] = N_toadd_01_J[k][filter]/N_01_J[k][filter]
			CF_10_J[k] = N_toadd_10_J[k][filter]/N_10_J[k][filter]
		endfor
		
		//Plot the correction factors
		Display/N=Fig6 CF_01_J
		Appendtograph CF_10_J
		ModifyGraph log(left)=1
		ModifyGraph mode=3,marker=19,rgb(CF_10_J)=(0,0,65280)
End


//ANCILLARY FUNCTIONS NOT CALLED BY THE END USER 
//(alphabetical order)

Function calcDudko([plot])
//Calculates the ideal curves based on the simulation parameters, k_01_ideal and k_10_ideal, 
//using the Dudko equation with nu=2/3.

	//Is a plot to be generated?
	Variable plot
	if(paramisDefault(plot))
		plot = 0
	endif
	
	//Reference to user-specified parameter wave
	Wave Param
	
	//Find the force range being probed
	Variable Z_start = param[4]
	Variable Z_stop = param[5]
	Variable ks = param[6]
	Variable Lc_linker = param[7]
	Variable kT = 4.1e-21
	Variable Lc_total = Lc_linker + param[8]
	Variable F_low = ks*(Z_start - calcXfromZ(Z_start, Lc_linker))
	Variable F_high = ks*(Z_stop - calcXfromZ(Z_stop, Lc_total)) 
	
	//Make the waves to store the output rate maps
	Make/O/D k_01_ideal, k_10_ideal
	SetScale/I x, F_low, F_high, k_01_ideal, k_10_ideal
	
	//Evaluate the Dudko equation:
	k_01_ideal = Param[11]*(1-2/3*x*Param[9]/Param[13])^(1/2)*exp(Param[13]/kT*(1-(1-2/3*x*Param[9]/Param[13])^(3/2)))
	k_10_ideal = Param[12]*(1-2/3*x*Param[10]/Param[14])^(1/2)*exp(Param[14]/kT*(1-(1-2/3*x*Param[10]/Param[14])^(3/2)))
	
	//If needed, generate an output plot:
	if(plot==1)
		Display k_01_ideal
		Appendtograph k_10_ideal
		ModifyGraph log(left)=1
	endif
End

Function calcXfromZ(Z, Lc)
//For particular polypeptide parameters and a particular cantilever position, determines the molecular extension.
//Z is the cantilever position and Lc is the polypeptide contour length.
	Variable Z, Lc
	
	//Load user-specified simulation parameters:
	Wave Param
	
	//Additional parameters
	Variable ks = Param[6]	//Spring constant
	Variable kT = 4.1e-21		//Thermal energy
	Variable lp = 0.4e-9		//Polypeptide persistence length
	
	//Solving for X as a numerical cubic-polynomial root-finding problem:
	Make/N=4/O polyWave
	polyWave[0] = -kT*Lc^2/(4*ks*lp)+Z*Lc^2
	polyWave[1] = -2*Z*Lc-Lc^2
	polyWave[2] = Z+2*Lc
	polyWave[3] = -1
	
	FindRoots/P=polyWave
	Wave W_polyRoots
	Variable X = real(W_polyRoots[0])
	
	//Cleanup:
	KillWaves polyWave, W_polyRoots
	
	//Return polymer extension:
	return X
End

Function estLoadingRate(force, Lc)
//Estimates the local loading rate (i.e. change in force per time) at a particular force for some contour length
//This is used in evaluating Eq. 3
	Variable force, Lc
	
	//Simulation parameters
	Wave Param
	
	//Parameters for this function
	Variable search_step = 1e-11
	Variable interval = 1e-10
	
	//Find the Z value (i.e., cantilever position) corresponding to the force in question
	Variable curr_Z = Param[4] //starting Z
	Variable SC = Param[6]
	Variable curr_F = 0
	do
		curr_Z += search_step
		curr_F = (curr_Z - calcXfromZ(curr_Z, Lc))*SC
	while (curr_F < force)
	Variable Z = curr_Z
	
	//Find the loading rate
	Variable interval_t = interval/Param[3]
	Variable low_F = ((Z-interval) - calcXfromZ(Z-interval, Lc))*SC
	Variable high_F = ((Z+interval) - calcXfromZ(Z+interval, Lc))*SC
	Variable loading_rate = (high_F - low_F)/(interval_t*2)
	
	//Return loading rate
	Return loading_rate
End

Function filterRecord(F_sim, state_sim, t_filter)
//Takes the simulation output (F_sim, state_sim) and applies a hard time-response filter, removing dwells 
//in states shorter than the t_filter
	Wave F_sim, state_sim
	Variable t_filter
	
	//Simulation parameters
	Wave Param
	
	//Convert t_filter from microseconds to seconds
	t_filter = t_filter*1e-6
	
	//Make waves to contain the filtered output
	Duplicate/O F_sim, F_sim_ds
	Duplicate/O state_sim, state_sim_ds
	
	//Initialize function parameters
	Variable curr_entry = 0
	Variable curr_state = state_sim[0]
	Variable dt = Param[0]
	
	//Iterate through each point in the simulated trace.  If a dwell is to be removed, change each force point
	//in the simulated record to the value corresponding to the other state
	Variable k, n, new_Lc, curr_Z, curr_X
	for (k=0; k<numpnts(F_sim_ds)-1; k+=1) //iteration over points in record
		if (state_sim_ds[k] != state_sim_ds[k+1]) //i.e., a transition
			if ( (k-curr_entry)*dt < t_filter) //dwell shorter than filter time and must be removed
				//Remove the dwell from the state wave:
				for (n=curr_entry; n<=k; n+=1)
					state_sim_ds[n] = state_sim_ds[k+1]
				endfor
				//Remove the dwell from the force wave (which involves a polymer elasticity calculation):
				if (state_sim_ds[k] == 1)
					new_Lc = Param[7] + Param[8]
				else
					new_Lc = Param[7]
				endif
				for (n=curr_entry; n<=k; n+=1)
					curr_Z = dimOffset(f_sim,0) + dimDelta(f_sim,0)*n
					curr_X = calcXfromZ(curr_Z, new_Lc)
					F_sim_ds[n] = (curr_Z-curr_X)*Param[6]	
				endfor
				
			else //dwell is ok
				curr_entry = k+1
			endif
			curr_state = state_sim_ds[k+1]
		endif
	endfor
End

Function simSingle(direction, [plot])
//Simulate a single force ramp as the AFM cantilever is moved either away from (direction = 1) or towards
//(direction = -1) the surface.  This function is repeatedly called and analyzed by simBatch() until the desired
//number of transitions are recorded.
	Variable direction, plot
	
	//Is a force-extension plot to be displayed
	if(paramisdefault(plot))
		plot = 0
	endif
	
	//Simulation parameters
	Wave param
	
	//Additional parameters parsed from the values in param
	Variable dt = param[0]
	Variable zstep = direction*param[3]*dt
	Variable num_steps = abs((param[5]-param[4])/zstep)
	Variable Lc_fold = param[7]
	Variable Lc_unfold = Lc_fold + param[8]
	Variable ks = param[6]
	Variable kT = 4.1e-21
	Variable dx_01 = Param[9]
	Variable k0_01 = Param[11]
	Variable dG_01 = Param[13]
	Variable dx_10 = Param[10]
	Variable k0_10 = Param[12]
	Variable dG_10 = Param[14]
	
	//Make waves that will store the output of the simulation:
	Make/O/N=(num_steps) F_sim, X_sim, state_sim
	
	//Initiate the simulation
	Variable Z, Lc
	if (direction == 1)
		Z = param[4]
		state_sim[0] = 0
		Lc = Lc_fold
	else 
		Z = param[5]
		state_sim[0] = 1
		Lc = Lc_unfold
	endif
	X_sim[0] = calcXfromZ(Z, Lc)
	F_sim[0] = ks*(Z-X_sim[0])
	
	SetScale/P x, Z, zstep, F_sim, X_sim, state_sim
	
	//Iterate through each time step in the simulation
	Variable k, curr_k, curr_P, curr_k0, curr_dx, curr_dG
	for (k=1; k<numpnts(X_sim); k+=1)
		//See if a transition is to occur:
		if (state_sim[k-1] == 0)
			curr_dx = dx_01
			curr_dG = dG_01
			curr_k0 = k0_01
		elseif (state_sim[k-1] == 1)
			curr_dx = dx_10
			curr_dG = dG_10
			curr_k0 = k0_10
		endif
		//Evaluation of rate at this force based on Dudko model:
		curr_k = curr_k0*(1-2/3*F_sim[k-1]*curr_dx/curr_dG)^(1/2)*exp(curr_dG/kT*(1-(1-2/3*F_sim[k-1]*curr_dx/curr_dG)^(3/2)))
		curr_P = curr_k*dt //Transition probability
		//Stochastic simulation to determine if a transition occurs:
		if (abs(enoise(1)) <= curr_P) //transition occurs
			state_sim[k] = 1-state_sim[k-1]
		else	//transition does not occur
			state_sim[k] = state_sim[k-1]
		endif
		
		//Increment values and compute new values
		Z += zstep
		if (state_sim[k] == 1)
			Lc = Lc_unfold
		else
			Lc = Lc_fold
		endif
		X_sim[k] = calcXfromZ(Z, Lc)
		F_sim[k] = ks*(Z-X_sim[k])
	
	endfor
	
	//Make a plot, if needed
	if (plot==1)
		display F_sim vs X_sim
		display F_sim
	endif
End