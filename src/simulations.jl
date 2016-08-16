#using ODE
#using Gillespie
using Sundials
using PyCall
#using RCall

@pyimport numpy.random as py_rand
#@rimport CircStats


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
1D diffusion gradient reaction function to integrate

Intput: 
	# x = cell position (as distance from wound)

	# t = absolute time after wounding (in seconds)

	# diffusion parameters:
	# D = diffusion coefficient
	# A = source strength
	# tau0 = time the wound is signalling (in seconds)
"""
function diffusion_gradient(tau::Float64,D::Float64,A::Float64,t::Float64,x::Float64)

	#x0 is the position of the wount - can be set differently if needed
	x0=0
	res = A*exp(-(x-x0)^2/(4*D*(t-tau)))/sqrt(4.0*D*(t-tau)*1pi)
	return res
end


"""
1D diffusion gradient reaction - solve the reaction's function
Intput: 
	# x = cell position (as distance from wound)
	# r = cell radios

	# t = absolute time after wounding (in seconds)

	# diffusion parameters:
	# D = diffusion coefficient
	# A = source strength
	# tau0 = time the wound is signalling (in seconds)
"""
function diffusion_gradient_solver(tau0::Float64,D::Float64,A::Float64,t::Float64,x::Float64,r::Float64)

	x1 = x-r
	x2 = x+r
	integrand1 = (tau) -> diffusion_gradient(tau,D,A,t,x1)
	integrand2 = (tau) -> diffusion_gradient(tau,D,A,t,x2)
	conc1 = quadgk(integrand1,0,min(tau0,t))
	conc2 = quadgk(integrand2,0,min(tau0,t))

	return (conc1[1],conc2[1])

end


"""
simulate a random walk for one cell - basic model of (w,b,p)
"""
function simulate_random_walk(n_steps::Int64,start_point::Tuple{Float64,Float64},w::Float64,p::Float64,b::Float64,bias_angle::Float64)

	#save two result vectors: step lengths and step angles
	steps_length = Array{Float64,1}(n_steps)
	steps_angles = Array{Float64,1}(n_steps)
	steps_points = Array{Tuple{Float64,Float64},1}(n_steps+1)
	steps_points[1]=start_point

	#define distributions

	#step length
	dt = 0.05
	truncated_normal = Truncated(Normal(0,1), 0, Inf)

	steps_length[1] , steps_angles[1] , steps_points[2] = simulate_random_walk_one_step(truncated_normal,dt,w,p,Inf,b,bias_angle,start_point)

	for i in 2:n_steps

		steps_length[i] , steps_angles[i] , steps_points[i+1] = simulate_random_walk_one_step(truncated_normal,dt,w,p,steps_angles[i-1],b,bias_angle,steps_points[i-1])	

	end

	return Random_Walk(steps_length,steps_angles,steps_points)

end

"""
simulate a random walk model where the bias is changed every time step
Input:
		- n_steps and b and bias_angle need to be of same length
"""
function simulate_random_walk_changing_bias(n_steps::Int64,start_point::Tuple{Float64,Float64},w::Float64,p::Float64,b::Array{Float64,1},bias_angle::Array{Float64,1})

	#save two result vectors: step lengths and step angles
	steps_length = Array{Float64,1}(n_steps)
	steps_angles = Array{Float64,1}(n_steps)
	steps_points = Array{Tuple{Float64,Float64},1}(n_steps+1)
	steps_points[1]=start_point

	#define distributions

	#step length
	dt = 0.001
	truncated_normal = Truncated(Normal(0,1), 0, Inf)

	steps_length[1] , steps_angles[1] , steps_points[2] = simulate_random_walk_one_step(truncated_normal,dt,w,p,Inf,b[1],bias_angle[1],start_point)

	for i in 2:n_steps
		steps_length[i] , steps_angles[i] , steps_points[i+1] = simulate_random_walk_one_step(truncated_normal,dt,w,p,steps_angles[i-1],b[i],bias_angle[i],steps_points[i-1])
	end

	return Random_Walk(steps_length,steps_angles,steps_points)

end

""" 
simulate one step of a random walk
"""
function simulate_random_walk_one_step(truncated_normal::Distributions.Truncated{Distributions.Normal{Int64},Distributions.Continuous},
										dt::Float64,
										w::Float64,
										p::Float64,
										persistant_angle::Float64,
										b::Float64,
										bias_angle::Float64,
										prev_point::Tuple{Float64,Float64})

	#randomly choose the next step length, according to the step distribution
	next_length = rand(truncated_normal) * sqrt(dt)

	#randomly choose if the next step is going to be biased or persistent
	rand_w = rand()

	if rand_w<=w
		#biased walk
		if b==0.0
			#sample angle from wrapped uniform distribution
			next_angle = pi*rand(Uniform(-1,1))
		elseif b==1.0
			next_angle = bias_angle
		else
			#sample angle from wrapped normal distribution with mean=bias direction and sigma=-2log(b)
			bias_kappa = 1/(-2*log(b))

			#generate next random angle by calling python's Von Mises distribution
			next_angle = py_rand.vonmises(bias_angle,bias_kappa)

			# """
			# #generate next random angle using Julia's distributions package
			# bias_d = VonMises(bias_angle,bias_kappa)
			# next_angle = rand(sampler(bias_d))
			# """
			
			# """
			# #generate the next random angle by calling R's wrapped normal distribution
			# next_angle = R"rwrpnorm(1, $bias_angle, $b)"[1]
			# #convert the angle (R returns values between 0 and 2pi)
			# next_angle = (next_angle + pi)%(2*pi)-pi
			# """
		end
	else
		#persistant walk
		if p==0.0
			#sample angle from wrapped uniform distribution
			next_angle = pi*rand(Uniform(-1,1))
		elseif p==1.0
			next_angle = persistant_angle
		else
			#sample angle from wrapped normal distribution with mean=prev_angle direction and sigma=-2log(b)
			if persistant_angle==Inf
				#there is no previous angle to base the distribution on
				next_angle = pi*rand(Uniform(-1,1))
			else
				persistant_kappa = 1/(-2*log(p))
				
				#generate next random angle by calling python's Von Mises distribution
				next_angle = py_rand.vonmises(persistant_angle,persistant_kappa)


				# """
				# #generate next random angle using Julia's distributions package
				# persistant_d = VonMises(steps_angles[i-1],persistant_kappa)
				# next_angle = rand(sampler(persistant_d))
				# """
				
				# """
				# #generate the next random angle by calling R's wrapped normal distribution
				# prev_angle = persistant_angle
				# next_angle = R"rwrpnorm(1, $prev_angle, $p)"[1]
				# #convert the angle (R returns values between 0 and 2pi)
				# next_angle = (next_angle + pi)%(2*pi)-pi
				# """
			end
		end
	end

	#TODO: check if needs to be - accodring to orientaiton along the x and y axis
	next_x = prev_point[1] + next_length*cos(next_angle)
	next_y = prev_point[2] + next_length*sin(next_angle)

	next_point = (next_x,next_y)

	return (next_length,next_angle,next_point)
end


"""
simulate a random walk model with a 1D diffusion gradient 
Input:
	n_steps : number of steps
	start_point : center of the cell 
	cell_r : cell radios
	bias_angle : the angle of the wound
	params: model parameters:
							- w - weight of bias (w-1 is weight of persistance)
							- mb - scaling parameter for bias
							- b0 - baseline observed bias
							- mp - scaling parameter for persistance
							- p0 - baseline observed persistance

							- A - concentration strength
							- D - diffusion coefficient
							- tau0 - time the wound is signalling

							- Kd - receptor binding parameter (dissociation constant)
							- R0 - number of receptors
"""
function simulate_random_walk_attractant_diffusion_gradient(n_steps::Int64,start_point::Tuple{Float64,Float64},cell_r::Float64,bias_angle::Float64,params::Array{Float64,1};start_time::Int64=0.0)

	#save two result vectors: step lengths and step angles
	steps_length = Array{Float64,1}(n_steps)
	steps_angles = Array{Float64,1}(n_steps)
	steps_points = Array{Tuple{Float64,Float64},1}(n_steps+1)
	steps_points[1]=start_point

	#save the ob and op values
	bias=Array{Float64,1}()
    persistance=Array{Float64,1}()

	#save parameters
	(w,mb,b0,mp,p0,A,D,tau0,Kd,R0) = params

	#step length distribution
	#TODO: fit a better one to the data
	dt = 0.001
	truncated_normal = Truncated(Normal(2,1), 0, Inf)
	#gamma_dist = Gamma(0.75,2.5) #control gamma fit
	#gamma_dist = Gamma(1.0,2.5) #DMSO gamma fit
	#gamma_dist = Gamma(1.0,3.0) #p38 gamma fit
	#gamma_dist = Gamma(0.7,3.2) #JNK gamma fit
	#gamma_dist = Gamma(0.9,3.5) #green data
	#gamma_dist = Gamma(0.65,3.5) #red data

	#each step is 30 seconds of time
	current_time = start_time

	for i in 1:n_steps

		current_time += 20.0 #each step is 20 seconds

		# if current_time>50400
		# 	println("time exceeded : ", current_time)
		# end

		#get next step size
		#next_length = rand(truncated_normal) * sqrt(dt)
		next_length = rand(truncated_normal)
		#next_length = rand(gamma_dist)

		steps_length[i] = next_length

		#randomly choose if the next step is going to be biased or persistent
		rand_w = rand()

		#the distance from the wound is the distance along the y position
		if i==1
			y=start_point[2]
		else
			y=steps_points[i-1][2]
		end

		#get concentration at front and back
		(af,ar) = diffusion_gradient_solver(tau0,D,A,current_time,y,cell_r)

		#get internal concentration of the cell
		Raf = receptor_attractant_conc(Kd,R0,af)
		Rar = receptor_attractant_conc(Kd,R0,ar)

		if (Raf-Rar)<0
		 	#DEBUG
		 	#println("NEGATIVE")
		end

		if rand_w<=w

			#observed bias
			ob = mb*(Raf-Rar)+b0
			#println(ob)
			push!(bias,ob)

			 if ob<0.0 || ob>1.0
			 	#the parameters should be rejected, so no walk will be generated
			 	#println("ob is ", ob, " -> discard params")
			 	return false
			 end

			#bias_angle = convert_y_angle_to_x(bias_angle)

			if ob==0.0
				#sample angle from wrapped uniform distribution
				next_angle = pi*rand(Uniform(-1,1))
			elseif ob==1.0
				#next angle is the bias angle
				next_angle = bias_angle
			else
				#get the next angle according to the observed bias and the von mises distribution
				bias_kappa = 1/(-2*log(ob))
				next_angle = py_rand.vonmises(bias_angle,bias_kappa)
			end

		else
			#observed persistant
			op = mp*(Raf-Rar)+p0
			push!(persistance,op)


			 if op<0.0 || op>1.0
			 	#the parameters should be rejected, so no walk will be generated
			 	#println("op is ", op, " -> discard params")
			 	return false
			 end


			if i==1
				prev_angle = pi*rand(Uniform(-1,1))
			else
				prev_angle = steps_angles[i-1]
			end

			#prev_angle = convert_y_angle_to_x(prev_angle)

			if op==0.0
				#sample angle from wrapped uniform distribution
				next_angle = pi*rand(Uniform(-1,1))
			elseif op==1.0
				#next angle is the bias angle
				next_angle = prev_angle
			else
				#get the next angle according to the observed bias and the von mises distribution
				persistant_kappa = 1/(-2*log(op))
				next_angle = py_rand.vonmises(prev_angle,persistant_kappa)
			end

		end

		#convert angle to the negative y-axis
		#TODO: check correctness!!!

		#steps_angles[i] = convert_x_angle_to_y(next_angle)
		steps_angles[i] = next_angle

		#get the next center position of the cell ( according to the angle with the x-axis)
		#TODO: check correctness!!!
		
		if i==1
			prev_point = start_point
		else
			prev_point = steps_points[i-1]
		end

		next_x = prev_point[1] + (next_length)*cos(next_angle)
		next_y = prev_point[2] + (next_length)*sin(next_angle)

		next_point = (next_x,next_y)
		steps_points[i+1] = next_point

		if next_y<0
			#simulation has reached the y-position of the wound, stop simulations
			return Random_Walk(steps_length[1:i],steps_angles[1:i],steps_points[1:i],bias,persistance)
		end


	end

	return Random_Walk(steps_length,steps_angles,steps_points,bias,persistance)

end


""" 
return the internal concentration of the cell
Input:
		- Kd - receptor binding parameter (dissociation constant)
		- R0 - number of receptors
		- ac - attractant concentration
"""
function receptor_attractant_conc(Kd::Float64,R0::Float64,ac::Float64)
	
	return 0.5*(Kd+R0+ac) - sqrt(0.25*((Kd+R0+ac)^2)-(R0*ac))
end

"""
convert the angle with the x-axis to an angle with the negative y-axis
"""
function convert_x_angle_to_y(x_angle::Float64)

	#TODO: FIX!

	y_angle = x_angle
	if y_angle>=0
		if y_angle<=pi/2
			y_angle += pi/2
		else
			y_angle = -x_angle
		end
	else
		if y_angle>=-pi/2
			y_angle = -x_angle
		else
			y_angle += pi/2
		end
	end

	return y_angle	

end

"""
convert an angle with the negative y-axis to the x-axis
"""
function convert_y_angle_to_x(y_angle::Float64)

	#TODO: FIX!

	x_angle = y_angle
	if x_angle>=0

		if x_angle<=pi/2
			x_angle = -y_angle
		else
			x_angle -= pi/2
		end
		
	else
		
		if x_angle>=-pi/2
			x_angle -= pi/2
		else
			x_angle = -y_angle
		end

	end

	return x_angle	

end
