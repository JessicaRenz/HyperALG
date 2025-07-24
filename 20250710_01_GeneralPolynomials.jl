using Oscar

cd(@__DIR__)  # setzt das Working Directory auf den Ordner des Skripts

function number2binary(n,L)
    binary = ""
    v = 2^(L-1)
    while v >= 1
        if n >= v 
            n = n-v 
            binary = binary * "1"
        else 
            binary = binary * "0"
        end
        v = v/2
    end
    while length(binary) < L 
        binary = "0" * binary 
    end
    return binary 
end

function binary2number(bin,L)
    num = 0 
    for i = 1:L 
        if bin[i] == '1' 
            num = num + 2^(L-i)
        end
    end
    return num
end

function possible_transitions(L)
    n_partners = []
    cumulative_partners = []
    partners = []
    for i = 0:2^L -1 
        vertex = number2binary(i,L)
        n_partners = push!(n_partners,0)
        for j = 1:L 
            if vertex[j] == '0'
                chars= collect(vertex) 
                chars[j] = '1'
                end_node = join(chars)
                end_node_int = binary2number(end_node,L)
                partners = push!(partners, end_node_int)
                n_partners[i+1] = n_partners[i+1] +1
            end
        end
        if i == 0
            cumulative_partners = push!(cumulative_partners,1)
        else 
            c = cumulative_partners[i] + n_partners[i] 
            #print("i: ", i, ", n_partners: ", n_partners[i],", cum_partners: ", cumulative_partners[i] )
            cumulative_partners = push!(cumulative_partners,c)
        end
    end
    return (n_partners = n_partners, cumulative_partners = cumulative_partners, partners = partners)
end

function edges(L,n_partners,cumulative_partners,partners)
    start=[]
    dest=[]
    s = 0
    j = 1
    while j<= n_partners[1]
        n = partners[j]
        start = push!(start,s)
        dest = push!(dest,n)
        j= j+1
    end
    for i = 1: 2^L -1
        j = 1
        while j <= n_partners[i+1]
            s = i 
            c = cumulative_partners[i+1]
            n = partners[c+j-1]
            start = push!(start,s)
            dest = push!(dest,n)
            j = j+1
        end
    end
    return(start = start, dest = dest)
end

function define_variables(L,start, dest)
    d = L*2^(L-1) -L 
    R, a, b = polynomial_ring(QQ, :a => 1:d, :b => 1:d)
    alist = []
    blist = []
    for i = 1:length(start)
        if start[i] != 0
            blist = push!(blist,i)
        end
        if dest[i] != 2^L-1
            alist = push!(alist,i)
        end
    end
    return (R, a, b , alist = alist, blist = blist)
end

function read_data(L,data_label)
    data = readlines(data_label)
    count = zeros(Int,2^L )
    for i = 1: length(data)
        bin = data[i]
        num = binary2number(bin,L)
        count[num + 1] = count[num+1]+1
    end
    n = length(data)
    N = count ./ n
    return N
end

function A_polynomials(L,R,N,start, dest, alist, blist)
    CFN = Vector{typeof(one(R))}(undef,2^L)
    CFN[1] = R(rationalize(N[1]))
    CBN = Vector{typeof(one(R))}(undef,2^L)
    CBN[2^L] = R(rationalize(N[2^L]))
    A = Vector{typeof(one(R))}(undef,2^L)
    A[1] = R(rationalize(1))
    A[2^L] = R(rationalize(1))
    for i = 1:(L-1)
        for n = 1:2^L-2
            bin = number2binary(n,L)
            num_1 = count(c -> c == '1',bin)
            if num_1 == i
                CFN[n+1] = R(rationalize(N[n+1]))
                for j = 1: length(dest)
                    if dest[j] == n 
                        s = start[j]
                        index = 0
                        for k = 1: length(alist)
                            if alist[k] == j 
                                index = k 
                            end
                        end
                        CFN[n+1]=CFN[n+1] + CFN[s+1]*a[index]
                    end
                end
            end
            if num_1 == L-i
                CBN[n+1] = R(rationalize(N[n+1]))
                for j = 1: length(start)
                    if start[j] == n 
                        g = dest[j]
                        index = 0
                        for k = 1: length(blist)
                            if blist[k]==j 
                                index = k 
                            end
                        end
                        CBN[n+1] = CBN[n+1] + CBN[g+1]*b[index]
                    end
                end
            end
        end
    end
    for n = 1: 2^L -2
        A[n+1] = CFN[n+1] + CBN[n+1] - R(rationalize(N[n+1]))
    end
    return A 
end

function polynomial_system(L,R,As,start, dest, alist, blist)
    P = Vector{typeof(one(R))}(undef,2^L-2)
    for n = 1: 2^L -2
        P[n] = -As[n+1]
        for i = 1:length(dest)
            if dest[i] == n 
                s = start[i]
                index = 0
                for k = 1: length(alist)
                    if alist[k] == i
                        index = k
                    end
                end 
                P[n] = P[n] + As[s+1]*a[index]
            end
        end
    end
    j = 0
    for n = 1:2^L -1
        bin = number2binary(n,L)
        num_1 = count(c -> c == '1',bin)
        if num_1 != 1
            sum = 0
            for i = 1: length(dest)
                if dest[i] == n 
                    s = start[i]
                    index = 0
                    for k = 1:length(alist)
                        if alist[k] == i 
                            index = k 
                        end
                    end
                    if index != 0
                        sum = sum + As[s+1]*a[index]
                    end
                end
            end
             for i = 1: length(dest)
                if dest[i] == n 
                    s = start[i]
                    index_a = 0
                    index_b = 0
                    for k = 1:length(alist)
                        if alist[k] == i 
                            index_a = k 
                        end
                    end
                    for l = 1:length(blist)
                        if blist[l] == i 
                            index_b = l 
                        end
                    end 
                    if index_a != 0 && index_b != 0
                        poly = - As[s+1]*a[index_a] + sum*b[index_b]
                        P = push!(P,poly)
                    end
                end
            end
        end
    end
    n = 2^L - 1
    sum = 0
    for i = 1:length(dest)
        if dest[i] == n 
            s = start[i] 
            sum = sum + As[s+1]
        end
    end
    for j = 1:length(dest)
        if dest[j] == n 
            s = start[j]
            index = 0
            for i = 1: length(blist)
                if blist[i] == j 
                    index = i 
                end
            end
            poly = -As[s+1] + sum*b[index]
            P = push!(P,poly)
        end
    end
    return P 
end

#Input data
L = 3
data_label = "toy.txt"

#-----main function starts----

transitions = possible_transitions(L)
n_partners = transitions.n_partners
cumulative_partners = transitions.cumulative_partners
partners = transitions.partners
edges_list = edges(L,n_partners,cumulative_partners,partners)
start = edges_list.start
dest = edges_list.dest
var = define_variables(L, start, dest)
N = read_data(L,data_label)
As = A_polynomials(L,var.R,N,start,dest,var.alist,var.blist)

P = polynomial_system(L,var.R,As,start,dest,var.alist,var.blist)
