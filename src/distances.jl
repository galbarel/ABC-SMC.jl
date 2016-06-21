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

	#get the distance for each array from the two datasets

	#println("size of d1 ", size(d1))
	#println("size of d2 ", size(d2))

	if ( size(d1)[1] != size(d2)[1] )
		error("euclidian_distance: data sets are of different dimentions")
		return
	end

	temp_matrix = Array{Array{Float64,1},1}(size(d1)[1])

	#println("d1 ", d1)
	#println("d2", d2)

	for i in 1:size(d1)[1]
		#get the distance between the vectors
		vec1 = d1[i]
		vec2 = d2[i]

		#println("vec1 ",vec1)
		#println("vec2 ", vec2)

		if ( size(vec1)[1] != size(vec2)[1] )
			error("euclidian_distance: one of the species has different dimentions in the two data sets")
			return
		end

		temp_matrix[i] = (vec1-vec2) .* (vec1-vec2)
	end

	#println(temp_matrix)

	return sqrt(sum(sum(temp_matrix)))
	
end