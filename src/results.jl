include("model.jl")
include("sampler.jl")
include("simulations.jl")
include("ABC_Gillespie.jl")
include("distances.jl")

using DataFrames
using Gadfly
using Compose
using Colors
using JLD
using MultivariateStats
using RCall
using EntropyEstimators


""" save results of ABC rejection or ABC SMC (each population is saved in a separated directory)

"""
function save_results(results::Array{ABC_population,1},model::ABC_Model,path::ASCIIString)

	#save the priors into a txt and jld file
    writedlm(string(path,"/priors.txt"),model.priors)
    jldopen(string(path,"/priors.jld"),"w") do file
        write(file,"priors",model.priors)
    end

    #save the particles from each population into a file
    for i in 1:length(results)

        #save_pop(i,results[i],model,path)
        compare_prior_posterior(path,pop_num=i)

    end
    
end

"""
save the results of one population 
"""
function save_pop(pop_num::Int64,pop::ABC_population,model::ABC_Model,path::ASCIIString)
    
    #create population directory (only if non existing)
    pop_dir = string(path,"population",pop_num)
    run(`mkdir -p $pop_dir`)

    #save the variables of the current population into a .jld file that can be later loaded
    jldopen(string(pop_dir,"/pop$pop_num.jld"),"w") do file
       write(file,"distances",pop.distances)
       write(file,"weights",pop.weights)
       write(file,"particles",pop.particles)
       write(file,"epsilon",pop.epsilon)
       write(file,"num_steps",pop.num_steps)
    end

    #open new file for writing
    pop_file = open(string(pop_dir,"/particles.txt"), "w")

    particles = pop.particles
    
    for j in 1:length(particles)
        write(pop_file,string(particles[j]), "\n")
    end

    close(pop_file)

    #save weights and particles
    weights_file = string(pop_dir,"/weights.txt")
    dis_file = string(pop_dir,"/distances.txt")
    writedlm(weights_file,pop.weights)
    writedlm(dis_file,pop.distances)
    

    df = DataFrame()
    all_plots = Array{Gadfly.Plot,1}(model.num_params*model.num_params)
    #for models with many parameters, create a plot with the histograms only
    hist_plots = Array{Gadfly.Plot,1}(model.num_params)
    df_cols = collect(1:model.num_params)

    hist_plot_num = 1
    scatter_plot_num = 1

    cols = round(Int,model.num_params/2)
    if (model.num_params%2)==1
        cols+=1
    end

    rows = 2

    if model.num_params==3
    	rows=1
    end

    for j in 1:model.num_params

        #save all the plots of the current parameter into one plot
        par_plots = Array{Gadfly.Plot,1}(model.num_params)

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
        all_plots[hist_plot_num] = plot_j
        hist_plots[j] = plot_j
        par_plots[j] = plot_j
        hist_plot_num += (model.num_params+1)

        #generate scatter plots with all other parameters
        for n in 1:model.num_params
            if n!=j
                
                #only create scatter plots with other parameters
                n_vals = [particles[k][n] for k in 1:length(particles)]
                scatter_j_n = plot(x=p_vals,y=n_vals,Geom.point,
                                   Guide.xlabel(string("Parameter ",j)),
                                   Guide.ylabel(string("Parameter ",n)))
                all_plots[scatter_plot_num] = scatter_j_n
                par_plots[n] = scatter_j_n
                scatter_plot_num +=1
            else
                scatter_plot_num +=1
            end
        end

        #save the plots of the current parameter into one plot
        cs = transpose(reshape([Context[render(par_plots[i]) for i in 1:length(par_plots)]], cols,rows))
        p = gridstack(cs)
        draw(PDF(string(pop_dir,"/parameter_",j,"_plots.pdf"),18inch,10inch),p)
        
    end

    #save the plots into one plot
    cs = transpose(reshape([Context[render(all_plots[i]) for i in 1:length(all_plots)]], model.num_params,model.num_params))
    p = gridstack(cs)
    draw(PDF(string(pop_dir,"/parameters_histograms.pdf"),12inch,12inch),p)

    #save only histogram plots
    cs = transpose(reshape([Context[render(hist_plots[i]) for i in 1:length(hist_plots)]], cols,rows))
    h = gridstack(cs)
    draw(PDF(string(pop_dir,"/parameters_histograms_only.pdf"),18inch,10inch),h)

    #save the dataframe into a file
    writetable(string(pop_dir,"/particles_df.csv"),df)

    #plot the weighted posteriror distributions
    plot_weighted_posterior(pop_dir)
    #TODO: change so that this is done during the loop, create a function for making the plots and call it twice

    #PCA analysis
    #robustness_analysis(pop_dir)
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

"""
plot random walks
"""
function plot_random_walks(walks::Array{Random_Walk,1},path::ASCIIString;walks_title="Random Walks")

    p_layers = Array{Layer,1}(length(walks))
    walk_color = distinguishable_colors(length(walks))

    for i in 1:length(walks)

        w = walks[i]
        points= w.steps_points

        point_x = [points[i][1] for i in 1:length(points)]
        point_y = [points[i][2] for i in 1:length(points)]

        #l = layer(x=point_x,y=point_y,Geom.line)
        l = layer(x=point_x,y=point_y,Geom.path(),Theme(default_color=walk_color[i]))
        p_layers[i] = l[1]
    end

    p = plot(p_layers,Guide.XLabel("X"),
                     Guide.YLabel("Y"),
                     Guide.Title(walks_title),
                     Coord.Cartesian(xmin=-5, xmax=5, ymin=-10, ymax=5)
                     )

                    #Coord.Cartesian(xmin=-5, xmax=5, ymin=-20, ymax=0)

    draw(PDF(string(path,"random_walks.pdf"),12inch,12inch),p)

    return p     

end

"""plotting function for transition matrix"""
function my_spy(M::AbstractMatrix, elements::Gadfly.ElementOrFunction...; mapping...)
    is, js, values = Gadfly._findnz(x->!isnan(x), M)
    n,m = size(M)
    df = DataFrames.DataFrame(i=is, j=js, value=values)
    plot(df, x=:j, y=:i, color=:value,
        Coord.cartesian(fixed=true, xmin=0.5, xmax=m+.5, ymin=0.5, ymax=n+.5),
        #Scale.color_continuous(minvalue=0.0, maxvalue=0.015),
        Scale.color_continuous,
        Geom.rectbin,
        Scale.x_continuous,
        Scale.y_continuous,
        elements...; mapping...)
end

"""
plot transition matrixes on the same figure
"""
function plot_TMs(TM_arr::Array{Array{Float64,2},1},path::ASCIIString,particle_num::Int64,pop_num::Int64;num_groups::Int64=4)

    #create population directory (only if non existing)
    pop_dir = string(path,"population",pop_num)
    run(`mkdir -p $pop_dir`)

    TM_plots = Array{Gadfly.Plot,1}(length(TM_arr))

    for i in 1:length(TM_arr)

        #plot the transition matrix
        T_plot = my_spy(TM_arr[i],Guide.xlabel(""),Guide.ylabel(""))
        TM_plots[i] = T_plot

    end

    #plot all the TMs in one plot
    cs = transpose(reshape([Context[render(TM_plots[i]) for i in 1:length(TM_plots)]], num_groups,num_groups))
    p = gridstack(cs)
    title = compose(context(0, 0, 1w, 0.25inch),text(0.5, 1.0, "Time after Injury / Distance from Wound", hcenter, vbottom))
    draw(PDF(string(pop_dir,"/TMs_particle",particle_num,".pdf"),12inch,12inch),vstack(title, p))

end


"""
robustness analysis - PCA on correlation matrix
"""
function robustness_analysis(path::ASCIIString)

    # particles = string(path,"particles_df.csv")

    # weights = string(path,"weights.txt")

    # #open the DF of the particles
    # post_dist_df = readtable(particles)

    # #open the DF of the weights
    # weights_df = readtable(weights,header=false)
    # weights_vec = Array(weights_df[:,1])

    # #get the wieghted posterior
    # post_dist_weighted = weighted_posterior(post_dist_df,weights_vec)
    

    # """
    # #get the correlation matrix from R - using the weights
    # cor_mat = rcopy(R"cov.wt($post_dist_df,wt=$weights_vec,cor=TRUE)$cor")
    
    # #run PCA using R
    # #cor_mat_pca_sdev = rcopy(R"princomp($cor_mat,cor=T)$sdev")
    # cor_mat_pca_sdev = rcopy(R"prcomp($cor_mat,cor=T)$sdev")

    # #cor_mat_pca_loadings = rcopy(R"princomp($cor_mat,cor=T)$loadings")
    # cor_mat_pca_loadings = rcopy(R"prcomp($cor_mat,cor=T)$rotation")

    # PCA_vars = cor_mat_pca_sdev

    # PCA_comps = cor_mat_pca_loadings
    # """

    #get the weighted posterior
    post_dist_weighted = readtable(string(path,"/weighted_particles_df.csv"))

    #get the correlation matrix
    post_dist_cor = cor(DataMatrix(post_dist_weighted))

    #run PCA on the correlation matrix
    post_dist_PCA = fit(PCA,Array(post_dist_cor))

    #get the variances of the components
    PCA_vars = principalvars(post_dist_PCA)

    #get the components
    PCA_comps = projection(post_dist_PCA)
   

    #plot the variances for the first 5 components
    comps = ["PC1","PC2","PC3","PC4","PC5"]
    p = plot(x=comps,y=PCA_vars[1:5],Geom.bar,
                     Guide.XLabel("Principal Components"),
                     Guide.YLabel("Variance"),
                     Guide.Title("Correlation Matrix PCA")
                     )

    draw(PDF(string(path,"/PCA_variance.pdf"),12inch,12inch),p)

    #parameters
    pars = ["w","mb","b0","mp","p0","A","D","tau0","Kd","R0"]

    #plot the first component for all parameters
    pc1=PCA_comps[:,1]
    p = plot(x=pars,y=pc1,Geom.bar, 
                     Guide.XLabel(""),
                     Guide.YLabel("PC1"),
                     Guide.Title("Sloppy")
                     )

    draw(PDF(string(path,"/PCA_PC1.pdf"),12inch,12inch),p)

    #plot the fifth component for all parameters
    pc5=PCA_comps[:,5]
    p = plot(x=pars,y=pc5,Geom.bar,
                     Guide.XLabel(""),
                     Guide.YLabel("PC5"),
                     Guide.Title("Stiff")
                     )

    draw(PDF(string(path,"/PCA_PC5.pdf"),12inch,12inch),p)

    """
    #get the projection
    amount = zeros(size(PCA_comps))
    for i in 1:size(PCA_comps)[1]
        #get the sum of the ith row
        #PC_sum = sum(PCA_comps[i,:].^2)
        for j in 1:size(PCA_comps)[2]
            #get the sum of the jth col
            PC_sum = sum(PCA_comps[:,j].^2)
            amount[i,j] = (PCA_comps[i,j]^2) / PC_sum
        end
    end

    #plot the first and fifth projections

    #plot the first component for all parameters
    pc1=amount[:,1]
    p = plot(x=pars,y=pc1,Geom.bar, 
                     Guide.XLabel(""),
                     Guide.YLabel("PC1"),
                     Guide.Title("Sloppy")
                     )

    draw(PDF(string(path,"PCA_PC1_projection.pdf"),12inch,12inch),p)

    #plot the fifth component for all parameters
    pc5=amount[:,5]
    p = plot(x=pars,y=pc5,Geom.bar,
                     Guide.XLabel(""),
                     Guide.YLabel("PC5"),
                     Guide.Title("Stiff")
                     )

    draw(PDF(string(path,"PCA_PC5_projection.pdf"),12inch,12inch),p)
    """

end


"""
sample from posterior population accodring to weights
"""
function weighted_posterior(particles::DataFrames.DataFrame,weights::Array{Float64,1})

        sampled_particles = DataFrame()
        
        for j in 1:size(particles)[2]
            #sample each of the parameters vector
            par_vec = particles[:,j]
            sampled = StatsBase.sample(par_vec,WeightVec(weights),size(particles)[1])
            insert!(sampled_particles,j,sampled,symbol("parameter$j"))
        end

        return sampled_particles
end

"""
plot weighted posterior distributions
"""
function plot_weighted_posterior(path::ASCIIString)

    #get the particles and weights
    particles = readtable(string(path,"/particles_df.csv"))
    weights_df = readtable(string(path,"/weights.txt"),header=false)
    weights= Array(weights_df[:,1])
    
    #get the weighted posterior
    weighted_particles = weighted_posterior(particles,weights)

    #save the weighted particles
    writetable(string(path,"/weighted_particles_df.csv"),weighted_particles)

    num_params = size(particles)[2]

    all_plots = Array{Gadfly.Plot,1}(num_params*num_params)
    #for models with many parameters, create a plot with the histograms only
    hist_plots = Array{Gadfly.Plot,1}(num_params)

    hist_plot_num = 1
    scatter_plot_num = 1

    cols = round(Int,num_params/2)
    if (num_params%2)==1
        cols+=1
    end

    rows=2
    if num_params==3
    	rows=1
    end

    for j in 1:num_params

        #save all the plots of the current parameter into one plot
        par_plots = Array{Gadfly.Plot,1}(num_params)

        #get all the current parameter values
        p_vals = weighted_particles[:,j]

        #generate a histogram plot
        plot_j = plot(x=p_vals,Geom.histogram(bincount=10),
                        Guide.xlabel(string("Parameter ",j)),
                        Guide.ylabel("Frequency"),
                        Guide.title("Weighted Parameters Histogram")
                        )
        all_plots[hist_plot_num] = plot_j
        hist_plots[j] = plot_j
        par_plots[j] = plot_j
        hist_plot_num += (num_params+1)

        #generate scatter plots with all other parameters
        for n in 1:num_params
            if n!=j
                
                #only create scatter plots with other parameters
                n_vals = weighted_particles[:,n]
                scatter_j_n = plot(x=p_vals,y=n_vals,Geom.point,
                                   Guide.xlabel(string("Parameter ",j)),
                                   Guide.ylabel(string("Parameter ",n)))
                all_plots[scatter_plot_num] = scatter_j_n
                par_plots[n] = scatter_j_n
                scatter_plot_num +=1
            else
                scatter_plot_num +=1
            end
        end

        #save the plots of the current parameter into one plot
        cs = transpose(reshape([Context[render(par_plots[i]) for i in 1:length(par_plots)]], cols,rows))
        p = gridstack(cs)
        draw(PDF(string(path,"/weighted_parameter_",j,"_plots.pdf"),18inch,10inch),p)
        
    end

    #save the plots into one plot
    cs = transpose(reshape([Context[render(all_plots[i]) for i in 1:length(all_plots)]], num_params,num_params))
    p = gridstack(cs)
    draw(PDF(string(path,"/weighted_parameters_histograms.pdf"),12inch,12inch),p)

    #save only histogram plots
    cs = transpose(reshape([Context[render(hist_plots[i]) for i in 1:length(hist_plots)]], cols,rows))
    h = gridstack(cs)
    draw(PDF(string(path,"/weighted_parameters_histograms_only.pdf"),18inch,10inch),h)

end


"""
compare prior and posterior distribution using Entropy
"""
function compare_prior_posterior(path::ASCIIString;pop_num=5,estimator="shrinkage")

    #get the prior distributions
    priors = jldopen(string(path,"/priors.jld"),"r") do file
        read(file, "priors")
    end

    pop_dir = string(path,"/population",pop_num,"/")

    #get the particles
    particles = readtable(string(pop_dir,"/particles_df.csv"))

    #get the weighted particles
    weighted_particles = readtable(string(pop_dir,"/weighted_particles_df.csv"))

    #TODO: check what happends with weighted particles

    entropy_arr = zeros(size(particles)[2])

    #go over each parameter
    for i in 1:size(particles)[2]
        #create the prior distribution
        prior_dist = rand(priors[i],size(particles)[1])

        #get the posterior distribution
        post_dist = particles[:,i]

        #get the entropy using the requsted estimate
        #prior_post_entropy = get_entropy(prior_dist, post_dist)
        prior_ent = get_entropy(prior_dist)
        post_ent = get_entropy(post_dist)
        entropy_arr[i] = prior_ent - post_ent

    end

    #plot the entropy values
    pars = ["w","mb","b0","mp","p0","A","D","tau0","Kd","R0"]
    p = plot(x=pars,y=entropy_arr,Geom.bar,
                     Guide.XLabel(""),
                     Guide.YLabel("Entropy"),
                     Guide.Title("Prior Posterior Entropy")
                     )

    draw(PDF(string(pop_dir,"/entropy.pdf"),12inch,12inch),p)

end