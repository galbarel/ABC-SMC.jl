using Distributions

"""
internal type of a random walk (for one cell)
"""
type Random_Walk

    steps_length::Array{Float64,1}
    steps_angles::Array{Float64,1}
    steps_points::Array{Tuple{Float64,Float64},1}

    #save bias and persistance parameters (relevant only for random_walk_gradient model)
    bias::Array{Float64,1}
    persistance::Array{Float64,1}

    """constructor for random walks without bias and persistance arrays"""
    function Random_Walk(steps_length::Array{Float64,1},steps_angles::Array{Float64,1},steps_points::Array{Tuple{Float64,Float64},1})
        new(steps_length,steps_angles,steps_points,zeros(1),zeros(1))
    end

    function Random_Walk(steps_length::Array{Float64,1},steps_angles::Array{Float64,1},steps_points::Array{Tuple{Float64,Float64},1},bias::Array{Float64,1},persistance::Array{Float64,1})
        new(steps_length,steps_angles,steps_points,bias,persistance)
    end

    function Random_Walk()
        new()
    end

end

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
    priors::Array{Union{Float64,Distributions.Uniform,Distributions.DiscreteUniform,Distributions.Normal,Distributions.LogNormal},1}

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
    priors::Array{Union{Float64,Distributions.Uniform,Distributions.DiscreteUniform,Distributions.Normal,Distributions.LogNormal},1}
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


type Random_Walk_model <: ABC_Model

    #model informaiton
    name::ASCIIString
    number::Int

    #parameters
    num_params::Int

    #prior distributions for parameters
    priors::Array{Union{Float64,Distributions.Uniform,Distributions.DiscreteUniform,Distributions.Normal,Distributions.LogNormal},1}
    #priors::Array{Any,1}

    #integration properties
    integration_mode::ASCIIString

    #walk properties - can be set so that there is only one type of walks, or so that there are different clusters of walks
    n_steps::Union{Int64,Array{Int64,1},Array{Array{Int64,1},1}} #number of steps for each walk - can all be the same, or each walk can have a different length
    n_walks::Int64 #number of total walks (cells)
    n_clusters::Int64 #number of spatio-temporal clusters
    start_point::Union{Array{Tuple{Float64,Float64},1},Array{Array{Tuple{Float64,Float64},1},1}}
    start_time::Union{Array{Int64,1},Array{Array{Int64,1},1}}
    cell_radius::Float64
    bias_angle::Float64

    """this is the model constructor"""
    function Random_Walk_model()

        new()

    end

end



""" internal representation of an Algorithm
"""
type ABC_Algorithm

    #data info
    time_points::FloatRange{Float64}
    data::Union{Array{Array{Float64,1},1},Array{Random_Walk,1},Array{Array{Random_Walk,1},1}}

    #epsilon
    epsilon::Array{Float64,1}

    #number of runs for each step
    particles::Int

    #perturbation kernel
    kernels::Array{Distributions.Uniform,1}

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
    num_steps::Int64
    epsilon::Float64

    #results
    data::Array{Array{Union{Array{Array{Float64,1},1},Array{Random_Walk,1}},1},1}
    particles::Array{Array{Float64,1},1}
    weights::Array{Float64,1}
    distances::Array{Float64,1}


end