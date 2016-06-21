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
    priors::Array{Any,1}
    #constant priors are saved as Float64
    #distribution priors are saved as one of three:
    # Distributions.Uniform
    # Distributions.Normal
    # Distributions.LogNormal
    # with the corresponding parameters

    #TODO: how to make array of specific distribution types
    #Array{Union{Distributions}}

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
    priors::Array{Any,1}
    #constant priors are saved as Float64
    #distribution priors are saved as one of three:
    # Distributions.Uniform
    # Distributions.Normal
    # Distributions.LogNormal
    # with the corresponding parameters

    #TODO: how to make array of specific distribution types
    #Array{Union{Distributions}}

    #integration properties
    integration_mode::ASCIIString
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

    #number of runs for each step
    particles::Int
    populations::Int

    ###TODO: check if number of populations is set according to epsilon

    #TODO: add kernel, number of simulations when using Gillespie (beta)

    """this is the model constructor"""
    function ABC_Algorithm()

        new()

    end

end