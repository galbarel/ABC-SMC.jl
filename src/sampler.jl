"""
ABC rejection for parameter inference

Algorithm:
	1. Sample parameters from prior distribtion
	2. simulate data with parameters accodring to given model
	3. if distance(data,simulated_data) < epsilon then accept the parameters
	4. repeat until reached the set number of repeats

Input:

	- alg: an internal type of the algorithm
	- model: an internal type of the model
	- distance_func: one of the possible: euclidian / 

Output:

	- posterior : Array of arrays, each one is the posterior distribtion for each parameter
	- counter : number of runs for the algorithm in total

"""
function abc_rejection(alg::ABC_Algorithm,
						model::ABC_Model,
						distance_func::ASCIIString)

	counter = 0 #general counter for the number of times the simulation runs
	num_accepted = 0 #counter for the number of accepted parameters
	posterior = Array{Array{Float64,1},1}(alg.particles) #output array, for each parameter we have a posterior distribution

	while num_accepted < alg.particles

		counter += 1 

		#sample a new parameter accodring to the prior
		params = sample_params_from_prior(model.priors)

		#simulate a dataset with the sampled parameters
		sim_data = run_simulation(alg,model,params)

		#check the distance between the original dataset and the simulated dataset
		distance = 0

		if distance_func == "euclidian"
			distance = euclidian_distance(alg.data,sim_data)
		else
			error("abc_rejection: the requested distance function: " ,distance_func, " is not implemented")
			return
		end

		if distance<=alg.epsilon[1] #there is only one epsilon
			#add to output vector
			num_accepted += 1
			posterior[num_accepted] = params
		end


	end

return (posterior,counter)

end

"""
ABC SMC for parameter inference
"""
function abc_smc()
end

"""
ABC SMC for model selection
"""
function abc_smc_model_selection()
end

"""
sample parameters from priors

Input:
	- priors (from the model) : an array with prior distribution (or constant Float64) for each parameter

Output:
	- params : array of Float64 with one sample per parameter
"""
function sample_params_from_prior(priors::Array{Any,1})

	params = zeros(length(priors))

	for i in 1:length(priors)

		if typeof(priors[i]) == Float64
			#constant prior
			params[i] = priors[i]

		else
			#prior is taken from the corresponding distribution
			params[i] = rand(priors[i],1)[1]

		end

	end

	return params

end

function run_simulation(alg::ABC_Algorithm,model::ABC_Model,params::Array{Float64,1})

	if model.integration_mode=="ODE"

		data = abc_ode(model.integration_function,model.initial_conditions,alg.time_points,params)

	elseif model.integration_mode=="Gillespie"

		#data = abc_gillespie(model.integration_function,stoch_matrix::Array{Int64,2},model.initial_conditions,time_end::Float64,params)

		data = Gillespie_abc(params,model.initial_conditions,model.integration_matrix,model.integration_hazards,alg.time_points)
	end

	return data

end