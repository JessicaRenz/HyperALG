using Oscar
using DelimitedFiles

cd(@__DIR__)  					# setting the Working Directory to the folder in which this script is stored

# converts a binary string into the corresponding number 
function binary2number(bin,L)
    num = 0 
    for i = 1:L 
        if bin[i] == '1' 
           num += 2^(L-i)
        end
    end
    return num
end

# figures out which nodes are directly connected
# the code structure of this function is taken from HyperHMM (https://github.com/StochasticBiology/hypercube-hmm/tree/main) and converted into Julia code
function possible_transitions(L)
    # n_partners[i] stores the number of states reachable from state i-1
    n_partners = zeros(Int, 2^L)
    
    # cumulative partners gives the starting index in 'partners' for the transitions corresponding to state i-1
    cumulative_partners = zeros(Int, 2^L)
    
    # partners stores all possible transitions for all states sequentially.
    # The transitions from state i-1 are located at indices cumulative_partners[i]:(cumulative_partners[i+1]-1) in 'partners'
    partners = Vector{Int}(undef,L*2^(L-1))			 
    idx = 1
    for i = 0:2^L -1 
       	n_partners[i+1] = 0
        for j = 1:L 
            if (i & (1 << (L-j))) == 0			# bit-wise comparison between the binary sequence for i and the sequence with only bit L-j is 1. 
            							#returns 1 only if in all comparison both bits are 1
                end_node_int = i | (1 << (L-j))		# flips bit L-j to a 1
                partners[idx] = end_node_int
                idx += 1
                n_partners[i+1] += 1
            end
        end
        if i == 0
            cumulative_partners[i+1] = 1
        else 
            c = cumulative_partners[i] + n_partners[i] 
            cumulative_partners[i+1] = c
        end
    end
    resize!(partners, idx-1)
    return (n_partners = n_partners, cumulative_partners = cumulative_partners, partners = partners)
end

# gives back a list of all edges
function edges(L,n_partners,cumulative_partners,partners)
    total_edges = sum(n_partners)
    start = Vector{Int}(undef, total_edges) 		# entry i stores the starting node of edge i
    dest = Vector{Int}(undef, total_edges)		# entry i stores the corresponding end node of edge i
    s = 0
    j = 1
    idx = 1
    for j in 1:n_partners[1] 				# finding all edges that start in node 0
        n = partners[j]
        start[idx] = s
        dest[idx] = n
        idx += 1
    end
    for i = 1: 2^L -1					# go through all remaining nodes and find the corresponding outgoing edges
        s=i
        c = cumulative_partners[i+1]
        for j in 1:n_partners[i+1]
            n = partners[c+j-1]
            start[idx] = s
            dest[idx] = n
            idx += 1
        end
    end
    return(start = start, dest = dest)
end

#defines the polynomial ring to work in 
function define_variables(L,start, dest)
    d = L*2^(L-1) -L 							# number of edges in the hypercube we need to assign variables (i.e. for which the weight is not always 1)
    R, a, b = polynomial_ring(QQ, :a => 1:d, :b => 1:d)			# defines an 'a' and a 'b' variable for every edge
    alist = Vector{Int}(undef, length(start))				# entry i indicates the number of the edge (=index in 'start' and 'dest') to which variable a[i] belongs
    blist = Vector{Int}(undef, length(start))				# entry i indicates the number of the edge (=index in 'start' and 'dest') to which variable b[i] belongs
    
    ai = 0
    bi = 0
    
    for i = 1:length(start)						# the edges starting in node 0 don't need a 'b' variable
        if start[i] != 0
            bi += 1
            blist[bi] = i
        end
        if dest[i] != 2^L-1						# the edges ending in the final node don't need an 'a' variable
            ai += 1
            alist[ai] = i
        end
    end
    
    resize!(alist, ai)
    resize!(blist, bi)
    
    return (R, a, b , alist = alist, blist = blist)
end

#read the data set in 
function read_data(L,data_label)
    N = zeros(Float64,2^L )						# entry i stores the proportion of node i-1 in the dataset
    n = 0
    for bin in eachline(data_label)
        num = binary2number(bin,L)
        N[num + 1] += 1
        n += 1
    end
    if n == 0
    	error("No data found in file $data_label")
    end
    N ./= n
    return N
end

# calculates the proportions of all trajectories that pass a certain node, taking into account forward and backward contributions
function A_polynomials(L,R,a,b,N,start, dest, alist, blist)
    #-------------------------------------------------------------------------------
    # Build adjacency lists from edge list (start, dest)
    # incoming[n+1] = indices of edge ending at node n
    # outgoing[n+1] = indices of edges starting from node n
    #-------------------------------------------------------------------------------
    incoming = [Int[] for _ in 1:2^L]
    outgoing = [Int[] for _ in 1:2^L]
    for j in 1:length(start)
    	push!(incoming[dest[j]+1],j)
    	push!(outgoing[start[j]+1],j)
    end
    
    #-----------------------------------------------------------------------------
    # Build lookup arrays for alist and blist
    # These map edge index j -> its position in alist/blist
    #-----------------------------------------------------------------------------
    a_index = zeros(Int, length(dest))
    for (k,j) in pairs(alist)
         a_index[j] = k
    end
    
    b_index = zeros(Int, length(start))
    for (k,j) in pairs(blist)
         b_index[j] = k
    end

    #------------------------------------------------------------------------------
    # Initialize coefficient arrays (as elements of the polynomial ring R)
    # CFN = "forward" coefficients
    # CBN = "backward" coefficients
    # A = final polynomial values (A[i] reports the proportion of trajectories that pass through node i-1)
    # Note: we use rationalize(N[...]) to convert floats to rationals before embedding them into the polynomial ring R
    #------------------------------------------------------------------------------
    CFN = fill(R(0),2^L)
    CFN[1] = R(rationalize(N[1]))
    
    CBN = fill(R(0),2^L)
    CBN[2^L] = R(rationalize(N[2^L]))
    A = fill(R(0),2^L)
    
    A[1] = R(rationalize(1))
    A[2^L] = R(rationalize(1))
    
    #------------------------------------------------------------------------------
    # Dynamic programming over the hypercube:
    # At step i, compute CFN for nodes with i ones in their binary label
    # At step L-i, compute CBN for nodes with L-i ones
    #------------------------------------------------------------------------------
    for i = 1:(L-1)
        for n = 1:2^L-2
            num_1 = count_ones(n)				# number of ones in a binary n
            
            #Forward recursion
            if num_1 == i
                CFN[n+1] = R(rationalize(N[n+1]))
                for j in incoming[n+1]				# edges ending at n
                    s = start[j]				# source node of edge j
                    CFN[n+1] += CFN[s+1]*a[a_index[j]]		# update forward coefficient
                end
            end
            
            # Backward recursion
            if num_1 == L-i
                CBN[n+1] = R(rationalize(N[n+1]))
                for j in outgoing[n+1]				# edges starting at n
                    g = dest[j]					# target node of edge j 
                    CBN[n+1] += CBN[g+1]*b[b_index[j]]		# update backward coefficient
                end
            end
        end
    end
    
    #-------------------------------------------------------------------------------
    # Combine forward and backward contributions into A
    #-------------------------------------------------------------------------------
    for n = 1: 2^L -2
        A[n+1] = CFN[n+1] + CBN[n+1] - R(rationalize(N[n+1]))
    end
    return A 
end

function polynomial_system(L,R,a,b,As,start, dest, alist, blist)
    #-------------------------------------------------------------------------------
    # build adjacency list
    #-------------------------------------------------------------------------------
    incoming = [Int[] for _ in 1:2^L] 				#edges ending at node n
    for j in 1:length(start)
    	push!(incoming[dest[j]+1],j)
    end
    
    #-------------------------------------------------------------------------------
    # Build lookup arrays (edge index -> position in alist/blist)
    #-------------------------------------------------------------------------------
    a_index = zeros(Int, length(dest))
    for (k,j) in pairs(alist)
    	a_index[j] = k
    end
    
    b_index = zeros(Int, length(start))
    for (k,j) in pairs(blist)
    	b_index[j] = k
    end
    
    #--------------------------------------------------------------------------------
    # Initialize polynomial list
    #--------------------------------------------------------------------------------
    P = Vector{typeof(one(R))}(undef,2^L-2)
    Pb_list =[]
    index_list = []
    
    #--------------------------------------------------------------------------------
    # First part: equations for the nodes from n=1,...,2^L -2
    #--------------------------------------------------------------------------------
    for n = 1: 2^L -2
        P[n] = -As[n+1]
        for j in incoming[n+1]			# edges ending at n 
            s = start[j]
            if a_index[j] != 0			# when a_index[j] == 0, this means the edge j has no corresponding value in alist
                P[n] += As[s+1]*a[a_index[j]]
            end
        end
    end
    
    #---------------------------------------------------------------------------------
    # Second part: constraints given through the introduction of the bi;j
    #---------------------------------------------------------------------------------
    for n = 1:2^L -1
        if count_ones(n) != 1
            sum = R(0)
            # sum over incoming edges 
            for j in incoming[n+1]		
                s = start[j]
                if a_index[j] !=0 
                     sum += As[s+1]*a[a_index[j]]
                end
            end
            
            # set the b_{i;j} into relation to the a_{i;j}
            for j in incoming[n+1]
                s = start[j]  
                if a_index[j] != 0 && b_index[j] != 0
                    poly = - As[s+1]*a[a_index[j]] + sum*b[b_index[j]]
                    push!(P,poly)
                    push!(Pb_list, length(P))
                    push!(index_list,b_index[j])
                end
            end
        end
    end
    
    #----------------------------------------------------------------------
    # Last node (all a_{i;n} are 1)
    #----------------------------------------------------------------------
    n = 2^L - 1
    sum = R(0)
    for j in incoming[n+1]
        s = start[j] 
        sum += As[s+1]
    end
    
    for j in incoming[n+1]
        s = start[j]
        if b_index[j] != 0
            poly = -As[s+1] + sum*b[b_index[j]]
            push!(P,poly)
            push!(Pb_list,length(P))
            push!(index_list,b_index[j])
        end
    end
    return (P_list =P,Pb_list =Pb_list, index_list =index_list) 
end

# removes one redundant polynomial for every step 1,...,L-1
function remove_polynomials(P,L)
    # store the last node index for each step in the hypercube
    last_node = zeros(Int,L-1)
    for j = 1:2^L -2
        w = count_ones(j)
        last_node[w] = j
    end
    
    # Build delete list 
    delete_list = filter(!=(0), last_node)
    deleteat!(P,delete_list)
    return P 
end

# eliminate redundant variables from the system
function remove_variables( start, dest, L, a ,alist, b, blist, P, Pb_list,index_list)
    #-------------------------------------------------------------------------------------------
    # Build adjacency lists
    # outgoing[n+1] = all edges starting at node n
    # incoming[n+1] = all edges ending at node n
    #--------------------------------------------------------------------------------------------
    outgoing = [Int[] for _ in 1:2^L]
    incoming = [Int[] for _ in 1:2^L]
    for e in 1:length(start)
    	push!(outgoing[start[e]+1], e)
    	push!(incoming[dest[e]+1], e)
    end
    
    #--------------------------------------------------------------------------------------------
    # Build lookup arrays: edge index -> position in alist/blist
    # If edge is not in alist/blist, entry stays 0
    #--------------------------------------------------------------------------------------------
    a_index = zeros(Int, length(dest))
    for (k,j) in pairs(alist)
    	a_index[j] = k
    end
    
    b_index = zeros(Int, length(start))
    for (k,j) in pairs(blist)
    	b_index[j] = k
    end
    
    #--------------------------------------------------------------------------------------------
    # Collect all variables into one vector [a;b]
    # length of a is stored so we can offset b correctly
    #--------------------------------------------------------------------------------------------
    vars_all = [a;b]
    length_a = length(a)
   
   #---------------------------------------------------------------------------------------------
   # Process a-variables:
   # For each node s, eliminate the last outgoing "a"-edge variable by rewriting it as 1 - sum(other outgoings a's)
   #---------------------------------------------------------------------------------------------
   for s in 0:2^L-2
   	edges = outgoing[s+1]
   	sum_a = 0
   	for (i,e) in enumerate(edges)
   		if a_index[e] == 0
   			continue
   		end
   		if i != length(edges)
   			sum_a += a[a_index[e]]
   		else 
   		 	vars_all[a_index[e]] = 1-sum_a
   		end
   	end
   end
   
   #-------------------------------------------------------------------------------------------
   # Process b-variables:
   # For each node n, eliminate the last incoming "b"-edge variable by rewriting it as 1 - sum(other incoming b's)
   # Keep track of which b-edges were eliminated  (b_red)
   #--------------------------------------------------------------------------------------------
    b_red = Int[]
    for n in 1:2^L-1
	edges = incoming[n+1]
	sum_b = 0
	for (i,e) in enumerate(edges)
		if b_index[e] == 0
			continue
		end
		if i != length(edges)
			sum_b += b[b_index[e]]
		else
			vars_all[length_a + b_index[e]] = 1- sum_b
			push!(b_red, e)
		end
	end
    end
  
  #----------------------------------------------------------------------------------------------
  # Evaluate all polynomials with the reduced variable set
  #-----------------------------------------------------------------------------------------------
    for k = 1: length(P)
        P[k] = evaluate(P[k],vars_all)
    end
    
    #---------------------------------------------------------------------------------------------
    # Remove polynomials corresponding to eliminated b-variables
    #---------------------------------------------------------------------------------------------
    delete_list = Int[]
    for e in b_red
        if b_index[e] == 0
        	continue
        end
        
        for k in 1:length(index_list)
		if index_list[k] == b_index[e]
                        push!(delete_list,Pb_list[k])
                        break
                end
        end
    end
    deleteat!(P,delete_list)
    return P 
end

# Get command line arguments
data_label = ARGS[1]

#Infer L from the first line of the file
first_line = open(data_label, "r") do io
	readline(io)
end

L = length(first_line)

println("Creating polynomials for L= $L with data from $data_label")

#-----main function starts----

transitions = possible_transitions(L);
n_partners = transitions.n_partners;
cumulative_partners = transitions.cumulative_partners;
partners = transitions.partners;
edges_list = edges(L,n_partners,cumulative_partners,partners);
start = edges_list.start;
dest = edges_list.dest;

# Give the memory of no longer needed variables free
transitions = nothing
n_partners = nothing
cumulative_partners = nothing
partners = nothing
edges_list = nothing

varis = define_variables(L, start, dest);
a = varis.a;
b = varis.b;
N = read_data(L,data_label);
As = A_polynomials(L,varis.R,a,b,N,start,dest,varis.alist,varis.blist);

N = nothing

P = polynomial_system(L,varis.R,a,b,As,start,dest,varis.alist,varis.blist);

As = nothing

P_red = remove_variables(start,dest,L,a,varis.alist,b,varis.blist,P.P_list,P.Pb_list,P.index_list);

P = nothing
P_final = remove_polynomials(P_red,L)

result_file = "polynomials_$data_label"
a_variables_file = "a_variables_$data_label"
b_variables_file = "b_variables_$data_label"

println("Writing polynomials in file $result_file and the variable - edges correspondences in $a_variables_file and $b_variables_file")

open(result_file, "w") do io
	for (i, val) in enumerate(P_final)
		println(io, "P[$i] = $val")
	end
end



header = [("variable", "from", "to")]

rows_a = [(a[i], start[varis.alist[i]], dest[varis.alist[i]]) for i in 1:length(varis.alist)]
writedlm(a_variables_file, rows_a)

rows_b = [(b[i], start[varis.blist[i]], dest[varis.blist[i]]) for i in 1:length(varis.blist)]
writedlm(b_variables_file, rows_b)
