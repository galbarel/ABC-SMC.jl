#using ODE
#using Gillespie
using Sundials
using PyCall

@pyimport numpy.random as py_rand

"""
deterministic simulations (ODE solver)
currently the function only uses the solver ode45 from the package ODE

input:
		- func: function to integrate
		- init_vals: an array of floats for the initial values of each species
		- time_points: a range of floats of the time points for the integration
		- params: an array of floats of the parameters for the integration

output:
		- an array of arrays, with the solutions for each species at each given time point
"""
function abc_ode(func::Function,init_vals::Array{Float64,1},time_points::FloatRange{Float64},params::Array{Float64,1})

	#solve using ODE package
	#convert to function with no params argument
	#g = (t,v) -> func(t,v,params)
	#(t,res) = ode45(g,init_vals,time_points,points=:specified,reltol = 1e-5, abstol = 1e-8)

	#solve using Sundails package
	g = (t,v,v_dot) -> func(t,v,v_dot,params)
	#println("values ",init_vals)
	#println("params ", params)
	d = Sundials.cvode(g,init_vals,collect(time_points))

	#convert Array{Float64,2} to Array{Array{Float64,1},1}
	res = Array{Array{Float64,1},1}(size(d)[1])
	for i in 1:size(d)[1]
		index=i
		arr = Array{Float64,1}(size(d)[2])
		for j in 1:size(d)[2]
			arr[j] = d[index]
	    	index+=size(d)[1]
	   	end
		res[i] = arr
    end

	return res

end

"""
stochastic simulations (SDE solver)
"""
function abc_sde()
	###TODO
end

"""
stochastic simulations (Gillespie Algorithm)

input:
		- func: function to integrate
		- stoch_matrix : ???
		- init_vals: an array of floats for the initial values of each species
		- time_end: last time point
		- params: an array of floats of the parameters for the integration

output:
		- an array of arrays, with the solutions for each species at each given time point

"""
function abc_gillespie(func::Function,stoch_matrix::Array{Int64,2},init_vals::Array{Int64,1},time_end::Float64,params::Array{Float64,1})

	result = ssa(init_vals, func, stoch_matrix, params, time_end)
	data = ssa_data(result)

	#TODO: convert data from  DataFrames.DataFrame to array of arrays so that the it includes only the species data points and is the same type as for the ODE return value

	return data

end

"""
simulate a random walk for one cell
"""
function simulate_random_walk(n_steps::Int64,start_point::Tuple{Float64,Float64},w::Float64,b::Float64,p::Float64,bias_angle::Float64)

	#save two result vectors: step lengths and step angles
	steps_length = Array{Float64,1}(n_steps)
	steps_angles = Array{Float64,1}(n_steps)
	steps_points = Array{Tuple{Float64,Float64},1}(n_steps)

	#define distributions

	#step length
	dt = 0.001
	truncated_normal = Truncated(Normal(0,1), 0, Inf)

	for i in 1:n_steps

		#randomly choose the next step length, according to the step distribution
		next_length = rand(truncated_normal) * sqrt(dt)
		steps_length[i] = next_length

		#randomly choose if the next step is going to be biased or persistent
		rand_w = rand()
		#next_angle = pi*rand(Uniform(-1,1))

		if rand_w<=w
			#biased walk
			if b==0.0
				#sample angle from wrapped uniform distribution
				next_angle = pi*rand(Uniform(-1,1))
			else
				#sample angle from wrapped normal distribution with mean=bias direction and sigma=-2log(b)
				bias_kappa = 1/sqrt(-2*log(b))
				"""bias_d = VonMises(bias_angle,bias_kappa)
				next_angle = rand(sampler(bias_d))"""
				next_angle = py_rand.vonmises(bias_angle,bias_kappa)
			end
		else
			#persistant walk
			if p==0.0
				#sample angle from wrapped uniform distribution
				next_angle = pi*rand(Uniform(-1,1))
			else
				#sample angle from wrapped normal distribution with mean=prev_angle direction and sigma=-2log(b)
				if i==1
					#there is no previous angle to base the distribution on
					next_angle = pi*rand(Uniform(-1,1))
				else
					persistant_kappa = 1/sqrt(-2*log(p))
					"""persistant_d = VonMises(steps_angles[i-1],persistant_kappa)
					next_angle = rand(sampler(persistant_d))"""
					next_angle = py_rand.vonmises(steps_angles[i-1],persistant_kappa)
				end
			end
		end

		steps_angles[i] = next_angle

		if i==1
			prev_point=start_point
		else
			prev_point=steps_points[i-1]
		end

		#TODO: check if needs to be - accodring to orientaiton along the x and y axis
		next_x = prev_point[1] + next_length*cos(next_angle)
		next_y = prev_point[2] + next_length*sin(next_angle)

		steps_points[i] = (next_x,next_y)

	end

	return Random_Walk(steps_length,steps_angles,steps_points)

end