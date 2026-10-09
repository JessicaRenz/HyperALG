include("src/HyperALG.jl")

cd(@__DIR__)  					# setting the Working Directory to the folder in which this script is stored

# Get command line arguments
input_path = ARGS[1]

data_label = basename(input_path)

#Infer L from the first line of the file
first_line = open(input_path, "r") do io
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
N = read_data(L,input_path);
As = A_polynomials(L,varis.R,a,b,N,start,dest,varis.alist,varis.blist);

N = nothing

P = polynomial_system(L,varis.R,a,b,As,start,dest,varis.alist,varis.blist);

As = nothing

P_red = remove_variables(start,dest,L,a,varis.alist,b,varis.blist,P.P_list,P.Pb_list,P.index_list);

P = nothing
P_final = remove_polynomials(P_red,L)

result_file = "output/polynomials_$data_label"
a_variables_file = "output/a_variables_$data_label"
b_variables_file = "output/b_variables_$data_label"

isempty("output/") || mkpath("output/")

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
