#using Distances

""" Euclidian distance between data sets
get the euclidian distance between 2 multidimentional arrays
both arrays need to have the same dimentions (and each array within each array should also have the same dimentions)

Input:
	- d1 - first array
	- d2 - second array

Output:
	- Distance in Float64 

"""
function euclidian_distance(d1::Array{Array{Float64,1},1},d2::Array{Array{Float64,1},1})

	#get the distance for each array from the two datasets)

	if ( size(d1)[1] != size(d2)[1] )
		error("euclidian_distance: data sets are of different dimentions")
		return
	end

	temp_matrix = Array{Array{Float64,1},1}(size(d1)[1])

	for i in 1:size(d1)[1]
		#get the distance between the vectors
		vec1 = d1[i]
		vec2 = d2[i]


		if ( size(vec1)[1] != size(vec2)[1] )
			error("euclidian_distance: one of the species has different dimentions in the two data sets")
			return
		end

		temp_matrix[i] = (vec1-vec2) .* (vec1-vec2)
	end

	return sqrt(sum(sum(temp_matrix)))
	
end

"""
Hellinger distance between 2 matrices (of same size)
"""
function hellinger_distance(m1::Array{Float64,2},m2::Array{Float64,2})
	
	d_sum = 0

	for i in 1:size(m1)[1]
		for j in 1:size(m1)[2]
			d_sum += (sqrt(m1[i,j])-sqrt(m2[i,j]))^2
		end
	end

	return (1/sqrt(2))*sqrt(d_sum)

end

"""
compute the transition matrix for an array of random walks (all need to have the same number of steps)
"""
function transition_matrix(walks::Array{Random_Walk,1};nrows::Int64=15,ncols::Int64=15)

	#create the matrix
	T = Array{Float64,2}(nrows,ncols)	

	nbins=15
	half_nbins=round(Int,nbins/2)

	#fill the matrix
	for i in 1:nrows
		
		"""
		#calculate i interval from 0 to 2pi and convert 
		#TODO: need to convert in the right order!!! from -pi to pi

		alpha_i_min , alpha_i_max = (2*pi*(i-1))/nbins , (2*pi*i)/nbins

		#convert to -pi to pi interval
		alpha_i_min = (alpha_i_min + pi)%(2*pi)-pi
		alpha_i_max = (alpha_i_max + pi)%(2*pi)-pi

		if i==15
			#TODO: find a way to correct in formula
			alpha_i_max = 1*pi
		end

		"""

		"""
		#calculate i interval from -pi to pi
		#TODO: find a way to make it 15 intervals
		if i<=half_nbins
			#the first half of the intervals goes from -pi to 0
			alpha_i_min , alpha_i_max = ((i-1-half_nbins)*pi)/half_nbins , ((i-half_nbins)*pi)/half_nbins

			if alpha_i_max==0.0
				#the middle interval goes from before 0.0 to after 0.0 (beacuse the number of intervals is 15)
				alpha_i_max = -alpha_i_min
			end
		else
			#the second half of the intervals goes from 0 to pi
			alpha_i_min , alpha_i_max = ((i-half_nbins)*pi)/half_nbins , ((i+1-half_nbins)*pi)/half_nbins
		end
		"""
		
		#16 intervals
		#alpha_i_min , alpha_i_max = ((i-1-half_nbins)*pi)/half_nbins , ((i-half_nbins)*pi)/half_nbins


		#calculate intervals from -pi to pi (15 intervals)
		if i-1==0
			alpha_i_min , alpha_i_max = -pi , -pi + (2pi/nbins)

		else
			alpha_i_min , alpha_i_max = -pi+(2pi/nbins*(i-1)) , -pi+(2pi/nbins*i)
		end

		#println("interval ", i ," ", alpha_i_min, ",", alpha_i_max)

		for j in 1:ncols

			"""
			#get the j interval
			alpha_j_min , alpha_j_max = (2*pi*(j-1))/nbins, (2*pi*j)/nbins
			#convert to -pi to pi interval
			alpha_j_min = (alpha_j_min + pi)%(2*pi)-pi
			alpha_j_max = (alpha_j_max + pi)%(2*pi)-pi

			if j==15
				#TODO: find a way to correct in formula
				alpha_j_max = 1*pi
			end
			"""

			"""
			if j<=half_nbins
				#the first half of the intervals goes from -pi to 0
				alpha_j_min , alpha_j_max = ((j-1-half_nbins)*pi)/half_nbins , ((j-half_nbins)*pi)/half_nbins

				if alpha_j_max==0.0
					#the middle interval goes from before 0.0 to after 0.0 (beacuse the number of intervals is 15)
					alpha_j_max = -alpha_j_min
				end
			else
				#the second half of the intervals goes from 0 to pi
				alpha_j_min , alpha_j_max = ((j-half_nbins)*pi)/half_nbins , ((j+1-half_nbins)*pi)/half_nbins
			end
			"""

			#16 intervals
			#alpha_j_min , alpha_j_max = ((j-1-half_nbins)*pi)/half_nbins , ((j-half_nbins)*pi)/half_nbins

			#calculate intervals from -pi to pi (15 intervals)
			if j-1==0
				alpha_j_min , alpha_j_max = -pi , -pi + (2pi/nbins)

			else
				alpha_j_min , alpha_j_max = -pi+(2pi/nbins*(j-1)) , -pi+(2pi/nbins*j)
			end

			#println("interval j: ", alpha_j_min, ",", alpha_j_max)


			current_sum = 0

			n_walks = length(walks)
			#n_steps = length(walks[1].steps_length)-1 #assuming each walk has the same number of steps
			all_steps=0

			for q in 1:n_walks
				#go over each walk
				angles = walks[q].steps_angles
				I_qr_ij = 0

				#in case number of stepts can change for each walk
				n_steps = length(angles)-1
				all_steps += n_steps

				for r in 1:n_steps #last time step is ignored because there is no next one for it
					#go over each step of the walk
					#get the angle for the current step and the next step
					alpha_q_r_i = angles[r+1]
					alpha_q_r_j = angles[r]

					if ( alpha_q_r_i>=alpha_i_min && alpha_q_r_i<alpha_i_max ) && ( alpha_q_r_j>=alpha_j_min &&  alpha_q_r_j<alpha_j_max )
						#both angles are within the interval
						I_qr_ij += 1
					end
				end #end r loop

				current_sum += I_qr_ij
			end #end q loop

			#current_sum = current_sum / (n_walks*n_steps) #assuming each walk has the same number of steps
			current_sum = current_sum / all_steps

			T[i,j] = current_sum

		end #end j loop
	end #end i loop

	return T

end