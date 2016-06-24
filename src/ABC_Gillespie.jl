using Distributions

"""
Gillespie algorithm for ABC SMC simulations

Input:
	- params : parameters values
	- init_values : inital values for the species
	- stoch_matrix : stochiometry matrix for the model
	- hazards : hazard functions of the model
	- times: time points for which simulated data is generated

Output:
	simulated data for all requested time points (in an array where each item is an array with data for all species)
"""
function Gillespie_abc(params::Array{Float64,1},init_values::Array{Float64,1},stoch_matrix::Array{Array{Int64,1},1},hazards::Array{Function,1},times::FloatRange{Float64})

	times = collect(times)
	time = 0.0
	time_counter = 1
	flag = true
	concentrations = Array{Array{Float64,1},1}(length(times))

	#TODO: make sure this is the way to deal with the wanted time points
	
	current_values = init_values

	#while simulation time hasn't exceeded the last time point
	while flag
		(next_values,next_time) = Gillespie_abc_one_step(params,current_values,stoch_matrix,hazards,time)
		concentrations[time_counter] = next_values
		time = next_time

		#as long as the next time point is larger than the next wanted time point, need to save the same concentration for all the time points in between
		while (time >=times[time_counter])
			time_counter +=1
			concentrations[time_counter] = next_values
			if time_counter >= length(times)
				flag = false
				break
			end
		end

	end

	return concentrations
end


"""
Gillespie algorithm for ABC SMC simulations - one step only

Input:
	- params : parameters values
	- values : current values for the species
	- stoch_matrix : stochiometry matrix for the model
	- hazards : hazard functions of the model
	- time : current time point

Output:
	- values : the next values for all species
	- t_next : next time point
"""
function Gillespie_abc_one_step(params::Array{Float64,1},values::Array{Float64,1},stoch_matrix::Array{Array{Int64,1},1},hazards::Array{Function,1},time::Float64)

	#get the hazard function for each reaction
	h_vals = zeros(length(hazards))
	for i in 1:length(hazards)
		#run the hazard function with its corresponding rate 
		h_vals[i] = hazards[i](params[i],values)
	end

	#get the combined hazard function
	h0 = sum(h_vals)

	#get the next time step
	#t_inc = pdf(Exponential(),h0)
	t_inc = (-1.0/h0) * log(rand(Uniform(0,1)))
	t_next = time + t_inc

	#get the probability for each reaction
	h_probs = h_vals ./ h0

	#choose a reaction according to their probabilities
	rand_int = rand(Uniform(0,1))
	cum_prob = 0.0
	chosen_hazard = 0

	for i in 1:length(hazards)

		chosen_hazard = i
		cum_prob = cum_prob + h_probs[i]
		if cum_prob > rand_int
			break
		end

	end

	#update values according to the chosen hazard
	#values is a one dimentional array of the current concentrations
	#we add to it the values of the array in the row that has the same index as the chosen hazard
	values = values .+ stoch_matrix[chosen_hazard]

	return (values,t_next)

end