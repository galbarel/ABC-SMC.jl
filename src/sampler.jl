using Distributions
using StatsBase

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
						model::ABC_Model)

	counter = 0 #general counter for the number of times the simulation runs
	num_accepted = 0 #counter for the number of accepted parameters
	posterior = Array{Array{Float64,1},1}(alg.particles) #output array, for each parameter we have a posterior distribution
	weights = Array{Float64,1}(alg.particles)

	while num_accepted < alg.particles

		counter += 1 

		#sample a new parameter accodring to the prior
		params = sample_params_from_prior(model.priors)

		#simulate a dataset with the sampled parameters and check if the parameter were accepted
		accepted_sim = simulate_and_compare_data(alg,model,params,1)

		if accepted_sim>0 #at least one simulated dataset was accepted
			#add to output vector
			num_accepted += 1
			posterior[num_accepted] = params
			#the weights in this case is the number of accepted simulations
			weights[num_accepted] = accepted_sim

			println("number accepted ", num_accepted)
		end

	end

	#normalize the weights
	weights = weights ./ sum(weights)


	return ABC_population(1,counter,alg.epsilon[1],posterior,weights)

end

"""
ABC SMC for parameter inference
"""
function abc_smc(alg::ABC_Algorithm,
					model::ABC_Model)

	#get epsilon values for each population
	epsilon_arr = alg.epsilon
	#TODO: add an option to run the simulation with automatically generated epsilons (in cases where the epsilon array in alg is empty)

	#save results from each population
	results = Array{ABC_population,1}(length(epsilon_arr))

	#run simulation for each population separately

	for pop in 1:length(epsilon_arr)

		println("simulating population: ", pop)

		if pop==1

			#the parameters will be sampled from the prior only
			#simply run ABC rejection once
			current_pop = abc_rejection(alg,model)

		else
			prev_pop = results[pop-1]
			current_pop = create_population(alg,model,pop,prev_pop.particles,prev_pop.weights)

		end

		#println(current_pop)
		results[pop] = current_pop

	end

	return results

end

"""
ABC SMC for model selection
"""
function abc_smc_model_selection()
end


########################################################## HELPER FUNCTIONS ######################################################################



"""
sample parameters from priors

Input:
	- priors (from the model) : an array with prior distribution (or constant Float64) for each parameter

Output:
	- params : array of Float64 with one sample per parameter
"""
function sample_params_from_prior(priors::Array{Union{Float64,Distributions.Uniform,Distributions.Normal,Distributions.LogNormal},1})

	params = zeros(length(priors))

	for i in 1:length(priors)

		if typeof(priors[i]) == Float64
			#constant prior
			params[i] = priors[i]

		else
			#prior is taken from the corresponding distribution
			params[i] = rand(priors[i])

		end

	end

	return params

end


"""
run simulation in order to generate simulated data set with the given parameters

Input:
	- alg : the internal type of the algorithm
	- model : the internal type of the model
	- params : the parameters to generate the dataset with

Output:

	- An array containing all generated datasets (one for deterministic, multiple for stochastic)
		* each elemnt in this array is an array of array, where each array includes data values for all species at a specific time time point

"""
function run_simulation(alg::ABC_Algorithm,model::ABC_Model,params::Array{Float64,1})

	#generate a number of data-sets (as set in the alg fields)
	if model.integration_mode=="ODE"
		integration_repeats=1
	else
		integration_repeats = model.integration_repeats
	end

	data = Array{Array{Array{Float64,1},1},1}(integration_repeats)

	if model.integration_mode=="ODE"

		#there is only one repeat for the integration
		res = abc_ode(model.integration_function,model.initial_conditions,alg.time_points,params)
		data[1] = res

	elseif model.integration_mode=="Gillespie"

		for i in 1:integration_repeats
			println("Gillespie run number: ", i)

			#data = abc_gillespie(model.integration_function,stoch_matrix::Array{Int64,2},model.initial_conditions,time_end::Float64,params)
			data_i = Gillespie_abc(params,model.initial_conditions,model.integration_matrix,model.integration_hazards,alg.time_points)
			data[i] = data_i
		end

	end

	return data

end


"""
create the next population for the ABC SMC method
"""
function create_population(alg::ABC_Algorithm,
							model::ABC_Model,
							pop_num::Int64,
							prev_particles::Array{Array{Float64,1},1},
							prev_weights::Array{Float64,1})

	particle = 0
	accepted_particles = 0
	particles = Array{Array{Float64,1},1}(alg.particles) #array of particles : each item is an array of parameters
	weights = Array{Float64,1}(alg.particles)

	while accepted_particles < alg.particles

		particle +=1

		#sample a new particle

		sampled_particles = sample_and_perturb(alg,model,prev_particles,prev_weights)
		#println("sampled particles ")
		#println(sampled_particles)

		#simulate data set with the sampled perturbed parameters
		num_accepted = simulate_and_compare_data(alg,model,sampled_particles,pop_num)

		if num_accepted>0 #at least one simulated dataset was accepted
			#add to output vector
			accepted_particles += 1
			particles[accepted_particles] = sampled_particles

			println("accepted particles ", accepted_particles)

			#calculate weight for the accepted particle
			p_weight = calculate_particle_weight(sampled_particles,model.priors,num_accepted,prev_weights,prev_particles,alg.kernels)
			weights[accepted_particles] = p_weight

		end

	end

	#normalise the weights
	weights = weights ./ sum(weights)

	#update the kernels according to the current population
	update_kernels(particles,alg,model)

	return ABC_population(pop_num,particle,alg.epsilon[pop_num],particles,weights)

end

"""
simulate a data set with the given parameters and comapre to the original dataset

Input:
	- alg : internal representation of the algorithm
	- model : internal representation of the model
	- params : the parameters for generating the simulated data set
"""
function simulate_and_compare_data(alg::ABC_Algorithm,
									model::ABC_Model,
									params::Array{Float64,1},
									pop_num::Int64)

	#simulate a dataset with the sampled parameters
	sim_data = run_simulation(alg,model,params)

	#save the distances for all simulation data
	distances = Array{Float64,1}(length(sim_data))

	#go over all data-sets (one for determinstic, multiple for stochastic)
	for i in length(sim_data)
		sim_data_i = sim_data[i]

		#check the distance between the original dataset and the current simulated dataset
		distance = 0

		if alg.distance_func == "euclidian"
			distance = euclidian_distance(alg.data,sim_data_i)
		else
			error("simulate_and_compare_data: the requested distance function: " ,alg.distance_func, " is not implemented")
			return
		end

		distances[i] = distance

	end

	#go over all the distances and count how many where accepted
	accepted_distances = 0

	for dis in distances

		if dis<=alg.epsilon[pop_num]
			accepted_distances +=1
		end

	end

	return accepted_distances
end



function calculate_particle_weight(particle::Array{Float64,1},
									priors::Array{Union{Float64,Distributions.Uniform,Distributions.Normal,Distributions.LogNormal},1},
									num_accepted::Int64,
									prev_weights::Array{Float64,1},
									prev_particles::Array{Array{Float64,1},1},
									kernels::Array{Distributions.Uniform,1})

	#get the prior of the particle - multiple all parameters priors
	p_prior = 1.0

	for i in 1:length(particle)

		param_prior = 1.0

		if typeof(priors[i]) == Float64
			#constant prior
				if particle[i]==priors[i]
					param_prior = 1.0
				else
					param_prior = 0.0
				end
		else
			param_prior = pdf(priors[i],particle[i])
		end

		p_prior = p_prior*param_prior
	end

	numerator = p_prior * num_accepted

	#go over all previous weights and multiple by the pertrubation kernel
	denominator = 0.0
	
	for j in 1:length(prev_particles)

		#go over each parameter and get the pdf of the kernel - in order to get the pdf of the particle
		kernel_prob = 1.0
		for k in 1:length(particle)

			#the parameters of the kernel distribution are now accodring to the previous parameter
			(k_min,k_max) = Distributions.params(kernels[k])
			current_kernel = Uniform(prev_particles[j][k] + k_min , prev_particles[j][k] + k_max)
			
			kernel_prob = kernel_prob * pdf(current_kernel,particle[k])
			

			"""kernel_prob = kernel_prob * pdf(kernels[k],particle[k])"""
		end

		denominator += (prev_weights[j] * kernel_prob)
	end
	
	return (numerator/denominator)

end


function sample_and_perturb(alg::ABC_Algorithm,
							model::ABC_Model,
							prev_particles::Array{Array{Float64,1},1},
							prev_weights::Array{Float64,1})

	zero_prior = true
	params = Array{Float64,1}()

	while zero_prior
		#as long as the prior for the perturb particle is zero, keeping sampling and perturbing

		"""
		#sample parameters from previous population 
		
		rand_u = rand(Uniform(0,1))
		cum_prob = 0

		#sum over all particles weights until reached the random number
		particle_index=0

		for p in 1:alg.particles
			#go over each particle
			#each particle is a vector of parameters
			#the entire vector is associated with a weight
			particle_index=p

			cum_prob = cum_prob + prev_weights[p]

			if (cum_prob>rand_u)
				break
			end
		end
		

		particle_index = p

		#particle_index is the chosen particle from the previous population

		sampled_from_prev = prev_particles[particle_index]
		"""

		sampled_from_prev = StatsBase.sample(prev_particles,WeightVec(prev_weights))

		#perturb the sampled parameters (accodring to kernel)
		#TODO: check correctness and add more kernels
		#the kernel is part of the alg and is a specified distribution
		#go over each parameter and perturb it (change in place of array)

		#prior probability of the particle
		particle_prior_prob = 1.0
		params = deepcopy(sampled_from_prev)

		#go over each parameter of the particle
		for i in 1:model.num_params
			
			params[i] += rand(alg.kernels[i])

			#get the prior probability of the perturbed parameter
			if typeof(model.priors[i]) == Float64
				#constant prior
				if params[i]==model.priors[i]
					param_prior_prob = 1.0
				else
					param_prior_prob = 0.0
				end
			else
				param_prior_prob = pdf(model.priors[i],params[i])
			end

			particle_prior_prob = particle_prior_prob * param_prior_prob

		end

		if particle_prior_prob>0
			zero_prior=false
		end
	end

	return params

end

"""
function to change (in place) perturbation kernels (a field in the algorim type) according to the previous population
"""
function update_kernels(particles::Array{Array{Float64,1},1},alg::ABC_Algorithm,model::ABC_Model)

	#go over each parameter in the population and update the kernel accordingly
	for i in 1:model.num_params
		#get all the current parameter values
		p_vals = [particles[j][i] for j in 1:length(particles)]
		scale = maximum(p_vals) - minimum(p_vals)
		alg.kernels[i] = Uniform(-scale/2.0,scale/2.0)
	end

end