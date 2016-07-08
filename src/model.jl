using Distributions
using DataFrames
using Gadfly
using Compose
using Colors

"""
internal type of a random walk (for one cell)
"""
type Random_Walk

    steps_length::Array{Float64,1}
    steps_angles::Array{Float64,1}
    steps_points::Array{Tuple{Float64,Float64},1}

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

    #walk properties
    n_steps::Int64 #number of steps for each walk
    n_walks::Int64 #number of total walks (cells)
    start_point::Tuple{Float64,Float64}
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
    data::Union{Array{Array{Float64,1},1},Array{Random_Walk,1}}

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
    particles::Array{Array{Float64,1},1}
    weights::Array{Float64,1}


end


""" save results of ABC rejection or ABC SMC (each population is saved in a separated directory)

"""
function save_results(results::Array{ABC_population,1},model::ABC_Model,path::ASCIIString)

    #save the particles from each population into a file
    for i in 1:length(results)

        #create population directory (only if non existing)
        pop_dir = string(path,"population",i)
        run(`mkdir -p $pop_dir`)

        #open new file for writing
        pop_file = open(string(pop_dir,"/particles.txt"), "w")

        pop = results[i]
        particles = pop.particles

        #go over each particle and add a line into the file
        #each row is a particle, each col is a parameter 
        #TODO: change so that each col is separated 
        
        for j in 1:length(particles)
            write(pop_file,string(particles[j]), "\n")
        end
        

        df = DataFrame()
        hist_plots = Array{Gadfly.Plot,1}(model.num_params*model.num_params)
        df_cols = collect(1:model.num_params)

        hist_plot_num = 1
        scatter_plot_num = 1

        for j in 1:model.num_params
            #get all the current parameter values
            p_vals = [particles[k][j] for k in 1:length(particles)]
            
            #add them to the dataframe
            insert!(df,j,p_vals,symbol("parameter$j"))

            #generate a histogram plot
            plot_j = plot(x=p_vals,Geom.histogram(bincount=10),
                            Guide.xlabel(string("Parameter ",j)),
                            Guide.ylabel("Frequency"),
                            Guide.title("Parameters Histogram")
                            )
            hist_plots[hist_plot_num] = plot_j
            hist_plot_num += (model.num_params+1)

            #generate scatter plots with all other parameters
            for n in 1:model.num_params
                if n!=j
                    
                    #only create scatter plots with other parameters
                    n_vals = [particles[k][n] for k in 1:length(particles)]
                    scatter_j_n = plot(x=p_vals,y=n_vals,Geom.point,
                                       Guide.xlabel(string("Parameter ",j)),
                                       Guide.ylabel(string("Parameter ",n)))
                    hist_plots[scatter_plot_num] = scatter_j_n
                    scatter_plot_num +=1
                else
                    scatter_plot_num +=1
                end
            end
            
        end

        #save the plots into one plot
        cs = transpose(reshape([Context[render(hist_plots[i]) for i in 1:length(hist_plots)]], model.num_params,model.num_params))
        p = gridstack(cs)
        draw(PDF(string(pop_dir,"/parameters_histograms.pdf"),12inch,12inch),p)

        #save the dataframe into a file
        writetable(string(pop_dir,"/particles_df.csv"),df)

        close(pop_file)

    end
    
end


"""
plot the data such that each variable is shown over time
"""
function plot_data(data::Array{Array{Float64,1},1},time_points::FloatRange{Float64},num_params::Int64,path::ASCIIString)
    
    df = DataFrame()
    df_cols = collect(1:num_params)
    p_layers = Layer[]
    p_colors = distinguishable_colors(num_params)
    times = collect(time_points)

    """p = plot(Guide.XLabel("Time"),
                 Guide.YLabel("Value"),
                 Guide.Title("Data"))
    """
    for j in 1:num_params
        #get all the current parameter values
        vals = [data[k][j] for k in 1:length(data)]
        
        #add them to the dataframe
        insert!(df,j,vals,symbol("var$j"))

        #create the layer
        l = layer(x=times,y=vals,Geom.point,Theme(default_point_size=8px,default_color=color(p_colors[j])))
        push!(p_layers,l[1])
        #Gadfly.add_plot_element!(p,l)
    end


    p = plot(p_layers,Guide.XLabel("Time"),
                     Guide.YLabel("Value"),
                     Guide.Title("Data"),
                     Guide.manual_color_key("Variable", [string("var",i) for i in 1:num_params], p_colors))

    
    draw(PDF(string(path,"data.pdf"),12inch,12inch),p)
    
end