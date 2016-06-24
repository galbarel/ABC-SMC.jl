using Distributions

abstract ABC_Model

""" internal representation of a model

#### Fields


#### Constructor
** Create a new model type **


"""
type Deterministic_Model <: ABC_Model

    #model informaiton
    name::ASCIIString
    number::Int

    #species and parameters
    num_species::Int
    num_params::Int

    #prior distributions for parameters
    priors::Array{Union{Float64,Distributions.Uniform,Distributions.Normal,Distributions.LogNormal},1}

    #integration properties
    integration_mode::ASCIIString
    integration_function::Function
    initial_conditions::Array{Float64,1}


    """this is the model constructor"""
    function Deterministic_Model()

        new()

    end

end

type Stochastic_Model <: ABC_Model

    #model informaiton
    name::ASCIIString
    number::Int

    #species and parameters
    num_species::Int
    num_params::Int

    #prior distributions for parameters
    priors::Array{Union{Float64,Distributions.Uniform,Distributions.Normal,Distributions.LogNormal},1}
    #priors::Array{Any,1}

    #integration properties
    integration_mode::ASCIIString
    integration_repeats::Int64
    integration_matrix::Array{Array{Int64,1},1}
    integration_hazards::Array{Function,1}
    initial_conditions::Array{Float64,1}

    """this is the model constructor"""
    function Stochastic_Model()

        new()

    end

end



""" internal representation of an Algorithm
"""
type ABC_Algorithm

    #data info
    time_points::FloatRange{Float64}
    data::Array{Array{Float64,1},1}

    #epsilon
    epsilon::Array{Float64,1}
    #TODO: set in the constractor - if the user does not input a list of epsilons, then generate an empty list which will indicate an automatic epsilon should be used

    #number of runs for each step
    particles::Int

    #perturbation kernel
    kernels::Array{Distributions.Uniform,1}
    #TODO: change type when more options for kernels are added

    #distance function
    distance_func::ASCIIString

    """this is the algorithm constructor"""
    function ABC_Algorithm()

        new()

    end

end

"""
internal type of a population from an ABC SMC simulation
"""
type ABC_population
    
    index::Int64

    #simulations information
    num_sampled::Int64
    epsilon::Float64

    #results
    particles::Array{Array{Float64,1},1}
    weights::Array{Float64,1}


end