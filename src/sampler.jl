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

Output:

	- internal type of ABC_population

"""
function abc_rejection(alg::ABC_Algorithm,
						model::ABC_Model,
						path::ASCIIString)

	counter = 0 #general counter for the number of times the simulation runs
	num_accepted = 0 #counter for the number of accepted parameters
	data = Array{Array{Union{Array{Array{Float64,1},1},Array{Random_Walk,1}},1},1}(alg.particles) #output array of the simulated data from each particle
	posterior = Array{Array{Float64,1},1}(alg.particles) #output array, for each parameter we have a posterior distribution
	weights = Array{Float64,1}(alg.particles)
	distances = Array{Float64,1}(alg.particles)

	while num_accepted < alg.particles

		counter += 1 

		#sample a new parameter accodring to the prior
		params = sample_params_from_prior(model.priors)

		#simulate a dataset with the sampled parameters and check if the parameter were accepted
		(accepted_sim,accepted_dis,sim_data,sim_TMs) = simulate_and_compare_data(alg,model,params,1)

		if sim_data==false
			#params need to be rejected
			continue
		end

		is_accepted = false
		#for the random walk model with gradient, each spatio-temporal cluster needs to have a distance < epsilon
		if model.integration_mode=="random_walk_gradient"
			"""
			if accepted_sim==model.n_clusters
				#accept params only is all clusters were below epsilon
				is_accepted=true
			end
			"""
			if sum(accepted_dis)<(model.n_clusters*alg.epsilon[1]) && accepted_sim>0
				#accept params when the sum of the distances was smaller than the sum of epsilons
				#and at least one simulation was accepted
				is_accepted=true

				#we want to plot the TMs 10% of the times
				if num_accepted%(alg.particles/10)==0 && num_accepted!=0
					plot_TMs(sim_TMs,path,num_accepted,1,num_groups=round(Int,sqrt(model.n_clusters)))
				end
			end
			
		else
			if accepted_sim>0
				#at least one simulated dataset was accepted
				is_accepted=true
			end
		end

		if is_accepted 
			#add to output vector
			num_accepted += 1
			data[num_accepted] = sim_data
			posterior[num_accepted] = params
			#the weights in this case is the number of accepted simulations
			weights[num_accepted] = accepted_sim
			distances[num_accepted] = maximum(accepted_dis)#accepted_dis[1] 
			#TODO: change that all distances are saved

			println("accepted particles: ", num_accepted)
			#println(" particles: ", params)
			#println("accepted distances: ", accepted_dis)

		end

	end

	#normalize the weights
	weights = weights ./ sum(weights)


	return [ABC_population(1,counter,alg.epsilon[1],data,posterior,weights,distances)]

end

"""
ABC SMC for parameter inference

Algorithm:
	1. create populations according to alg specifications
	2. for the first population : run ABC rejection algorithm with the first epsilon
	3. for any other population:
		3.1. sample particle from the previous population
		3.2 perturb particle according to a pertrubation kernel
		3.3. simulate a data set using the perturbed particle and compare with original data set with the corresponding epsilon
		3.4 if particle is accepted, calculate its weight
	4. return all populations

Input:

	- alg: an internal type of the algorithm
	- model: an internal type of the model
	- pop_num : number of populations when using adaptive epsilons

Output:

	- internal type of ABC_population
"""
function abc_smc(alg::ABC_Algorithm,
					model::ABC_Model,
					path::ASCIIString;
					pop_num::Int64=5)

	#get epsilon values for each population
	epsilon_arr = alg.epsilon

	num_pop = length(alg.epsilon)

	#adaptive epsilon
	if length(alg.epsilon)==1
		num_pop = pop_num
	end

	#save results from each population
	results = Array{ABC_population,1}(num_pop)

	#run simulation for each population separately

	for pop in 1:num_pop

		println("population number: ", pop)
		tic()

		if pop==1

			#the parameters will be sampled from the prior only
			#simply run ABC rejection once
			current_pop = abc_rejection(alg,model,path)[1]

		else
			prev_pop = results[pop-1]
			current_pop = create_population(alg,model,pop,prev_pop.particles,prev_pop.weights,prev_pop.distances,path)

		end

		results[pop] = current_pop
		#save the results of the current population
		save_pop(pop,current_pop,model,path)
		toc()

	end

	return results

end

"""
ABC SMC for model selection
"""
function abc_smc_model_selection()
#TODO
end


########################################################## HELPER FUNCTIONS ######################################################################



"""
sample parameters from priors

Input:
	- priors (from the model) : an array with prior distribution (or constant Float64) for each parameter

Output:
	- params : array of Float64 with one sample per parameter
"""
function sample_params_from_prior(priors::Array{Union{Float64,Distributions.Uniform,Distributions.DiscreteUniform,Distributions.Normal,Distributions.LogNormal},1})

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
	if model.integration_mode=="Gillespie"
		integration_repeats = model.integration_repeats
	elseif model.integration_mode=="random_walk_gradient"
		integration_repeats= model.n_clusters
	else
		integration_repeats=1
	end

	data = Array{Union{Array{Array{Float64,1},1},Array{Random_Walk,1}},1}(integration_repeats)

	if model.integration_mode=="ODE"

		#there is only one repeat for the integration
		res = abc_ode(model.integration_function,model.initial_conditions,alg.time_points,params)
		data[1] = res

	elseif model.integration_mode=="Gillespie"

		for i in 1:integration_repeats

			#println("Gillespie run ", i)

			#data_i = abc_gillespie(model.integration_function,stoch_matrix::Array{Int64,2},model.initial_conditions,time_end::Float64,params)
			#println("params : ", params)
			data_i = Gillespie_abc(params,model.initial_conditions,model.integration_matrix,model.integration_hazards,alg.time_points)
			data[i] = data_i
		end

	elseif model.integration_mode=="random_walk"

		walks = Array{Random_Walk,1}(model.n_walks)

		for i in 1:model.n_walks
			#simulating random walks (BM)
			walks[i] = simulate_random_walk(model.n_steps,model.start_point,params[1],params[2],params[3],model.bias_angle)
		end

		data[1]=walks

	elseif model.integration_mode=="random_walk_gradient"

		for j in 1:model.n_clusters

			walks = Array{Random_Walk,1}(model.n_walks)

			for i in 1:model.n_walks
				#simulating random walks (BM)
				current_walk = simulate_random_walk_attractant_diffusion_gradient(model.n_steps[j][i],model.start_point[j][i],model.cell_radius,model.bias_angle,params,start_time=model.start_time[j][i])
				if current_walk==false
					#could not run the simulation with the given parameters -> parameters should be discarded
					return false
				else
					walks[i] = current_walk
				end
			end

			data[j]=walks

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
							prev_weights::Array{Float64,1},
							prev_dis::Array{Float64,1},
							path::ASCIIString)

	particle = 0
	accepted_particles = 0
	data = Array{Array{Union{Array{Array{Float64,1},1},Array{Random_Walk,1}},1},1}(alg.particles) #output array of the simulated data from each particle
	particles = Array{Array{Float64,1},1}(alg.particles) #array of particles : each item is an array of parameters
	weights = Array{Float64,1}(alg.particles)
	distances = Array{Float64,1}(alg.particles)

	#check if epsilon needs to be updated according to the previous population
	if length(alg.epsilon)<pop_num
		#in this case the epsilon for populations >1 is taken to be the 10% quartile of the distances from the previous population
		alg.epsilon = push!(alg.epsilon,percentile(prev_dis,10))
		println("next epsilon ", alg.epsilon[pop_num])
	end

	while accepted_particles < alg.particles

		particle +=1

		#sample a new particle

		sampled_particles = sample_and_perturb(alg,model,prev_particles,prev_weights)

		#simulate data set with the sampled perturbed parameters
		(num_accepted,accepted_dis,sim_data,sim_TMs) = simulate_and_compare_data(alg,model,sampled_particles,pop_num)

		if sim_data==false
			#params need to be rejected
			continue
		end

		is_accepted = false
		#for the random walk model with gradient, each spatio-temporal cluster needs to have a distance < epsilon
		if model.integration_mode=="random_walk_gradient"
			"""
			if num_accepted==model.n_clusters
				#all clusters were accepted
				is_accepted=true
			end
			"""
			if sum(accepted_dis)<(model.n_clusters*alg.epsilon[pop_num]) && num_accepted>0
				#accept params when the sum of the distances was smaller than the sum of epsilons
				is_accepted=true

				#we want to plot the TMs 10% of the times
				if accepted_particles%(alg.particles/10)==0 && accepted_particles!=0
					plot_TMs(sim_TMs,path,accepted_particles,pop_num,num_groups=round(Int,sqrt(model.n_clusters)))
				end
			end
			
		else
			if num_accepted>0
				#at least one simulated dataset was accepted
				is_accepted=true
			end
		end

		if is_accepted #at least one simulated dataset was accepted
			#add to output vector
			accepted_particles += 1
			data[accepted_particles] = sim_data
			particles[accepted_particles] = sampled_particles

			#calculate weight for the accepted particle
			p_weight = calculate_particle_weight(sampled_particles,model.priors,num_accepted,prev_weights,prev_particles,alg.kernels)
			weights[accepted_particles] = p_weight
			distances[accepted_particles] = maximum(accepted_dis)#accepted_dis[1]

			println("accepted particles: ", accepted_particles)
			#println(" particles: ", sampled_particles)

		end

	end

	#normalise the weights
	weights = weights ./ sum(weights)

	#update the kernels according to the current population
	update_kernels(particles,alg,model)

	return ABC_population(pop_num,particle,alg.epsilon[pop_num],data,particles,weights,distances)

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

	if sim_data==false
		#parameters needs to be discarded
		return (0,0,sim_data,0)
	end

	#save the distances for all simulation data
	distances = Array{Float64,1}(length(sim_data))

	#save the transition matrix of the simulated data
	simulated_TMs = Array{Array{Float64,2},1}(length(sim_data))

	#go over all data-sets (one for determinstic, multiple for stochastic)
	for i in 1:length(sim_data)

		sim_data_i = sim_data[i]

		#check the distance between the original dataset and the current simulated dataset
		distance = 0

		if alg.distance_func == "euclidian"
			distance = euclidian_distance(alg.data,sim_data_i)
		elseif alg.distance_func == "hellinger"
			#check the distance between each spatio-temporal cluster separately
			sim_data_TM = transition_matrix(sim_data_i)
			simulated_TMs[i] = sim_data_TM
			distance = hellinger_distance(transition_matrix(alg.data[i]),sim_data_TM)
		else
			error("simulate_and_compare_data: the requested distance function: " ,alg.distance_func, " is not implemented")
			return
		end

		distances[i] = distance

	end

	#go over all the distances and count how many where accepted
	num_accepted_distances = 0
	#accepted_dis = 0
	accepted_dis = Array{Float64,1}()

	for dis in distances

		if dis<=alg.epsilon[pop_num]
			num_accepted_distances +=1
			#accepted_dis=dis
			push!(accepted_dis,dis)
		end

	end

	return (num_accepted_distances,accepted_dis,sim_data,simulated_TMs)
end



function calculate_particle_weight(particle::Array{Float64,1},
									priors::Array{Union{Float64,Distributions.Uniform,Distributions.DiscreteUniform,Distributions.Normal,Distributions.LogNormal},1},
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
		#the kernel is part of the alg and is a specified distribution
		#go over each parameter and perturb it (change in place of array)

		#prior probability of the particle
		particle_prior_prob = 1.0
		params = deepcopy(sampled_from_prev)

		#go over each parameter of the particle
		for i in 1:model.num_params
			
			params[i] += rand(alg.kernels[i])

			if typeof(model.priors[i]) == Distributions.DiscreteUniform
				#need to round the value of the perturbed parameter
				params[i] = round(Int,params[i])
			end

			#get the prior probability of the perturbed parameter
			if typeof(model.priors[i]) == Float64
				#round the value
				params[i] = round(params[i])
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
