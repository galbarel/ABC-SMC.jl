using ODE
#using Gillespie
#using Sundials

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

	#convert to function with no params argument
	g = (t,v) -> func(t,v,params)

	(t,res) = ode45(g,init_vals,time_points,points=:specified)

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