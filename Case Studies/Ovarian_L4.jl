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

function A_polynomials(L,R,a,b,N,start, dest, alist, blist)
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
                                break
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
                                break
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

function polynomial_system(L,R,a,b,As,start, dest, alist, blist)
    P = Vector{typeof(one(R))}(undef,2^L-2)
    Pb_list =[]
    index_list = []
    for n = 1: 2^L -2
        P[n] = -As[n+1]
        for i = 1:length(dest)
            if dest[i] == n 
                s = start[i]
                index = 0
                for k = 1: length(alist)
                    if alist[k] == i
                        index = k
                        break
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
                            break
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
                            break
                        end
                    end
                    for l = 1:length(blist)
                        if blist[l] == i 
                            index_b = l 
                            break
                        end
                    end 
                    if index_a != 0 && index_b != 0
                        poly = - As[s+1]*a[index_a] + sum*b[index_b]
                        #println("start: $s , end: $n , index_b: $index_b, poly: $poly")
                        P = push!(P,poly)
                        length_P = length(P)
                        #println("Index in P: $length_P")
                        Pb_list= push!(Pb_list, length_P)
                        index_list = push!(index_list,index_b)
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
                    break
                end
            end
            poly = -As[s+1] + sum*b[index]
            P = push!(P,poly)
            length_P = length(P)
            Pb_list = push!(Pb_list,length_P)
            index_list = push!(index_list,index)
        end
    end
    return (P_list =P,Pb_list =Pb_list, index_list =index_list) 
end

function remove_polynomials(P,L)
    delete_list = []
    for i = 1:L
        number = 0
        for j = 1:2^L -2
            bin = number2binary(j,L)
            num_1 = count(c -> c == '1',bin) 
            if num_1 == i
                number = number +1
            end
        end
        count_nodes = 0
        for j = 1:2^L-2
            bin = number2binary(j,L)
            num_1 = count(c -> c == '1',bin)
            if num_1 == i 
                count_nodes = count_nodes +1
                if count_nodes == number 
                    delete_list = push!(delete_list,j)
                end 
            end
        end
    end
    println("delete_list_2: $delete_list")
    deleteat!(P,delete_list)
    return P 
end

function remove_variables(n_partners, start, dest, L, a ,alist, b, blist, P, Pb_list,index_list)
    nb = zeros(Int,2^L)
    na = zeros(Int,2^L)
    for n = 0:2^L-1
        na[n+1] = n_partners[n+1]
        cnt = 0
        for i = 1:length(dest)
            if dest[i] == n 
                cnt = cnt +1
            end
        end
        nb[n+1] = cnt
    end
    vars_all = [a;b]
    length_a = length(a)
    for s = 0:2^L-2
        i = 0
        sum_a = 0
        for e = 1:length(start)
            if start[e] == s 
                i = i+1
                if i != na[s+1]
                    add_a = 0
                    for j = 1:length(alist)
                        if alist[j] == e 
                            add_a = add_a + a[j]
                            break
                        end
                    end
                    sum_a = sum_a + add_a 
                else
                    for j = 1:length(alist)
                        if alist[j] == e 
                            vars_all[j] = 1 - sum_a
                            break
                        end
                    end
                end
            end
        end
    end
    b_red = []
    for n = 1:2^L-1
        i = 0
        sum_b = 0
        for e = 1:length(dest)
            if dest[e] == n 
                i = i+1
                if i != nb[n+1]
                    add_b = 0
                    for j = 1:length(blist)
                        if blist[j] == e 
                            add_b = add_b + b[j]
                            break
                        end
                    end
                    sum_b = sum_b + add_b 
                else 
                    for j = 1:length(blist)
                        if blist[j] == e 
                            vars_all[length_a+j] = 1 - sum_b
                            append!(b_red,e)
                            break
                        end
                    end
                end
            end
        end
    end
    println("vars_all: $vars_all")
    for k = 1: length(P)
        P[k] = evaluate(P[k],vars_all)
    end
    println("P[29]: ")
    print(P[29])
    delete_list=[]
    for i = 1:length(b_red)
        P_entry = 0
        e = b_red[i]
        for j = 1:length(blist)
            if blist[j] == e 
                b_index = j
                for k = 1:length(index_list)
                    if index_list[k] == b_index
                        P_entry = Pb_list[k]
                        delete_list = push!(delete_list,P_entry)
                        println("delete_list: $delete_list")
                        break
                    end
                end
                break
            end
        end
    end
    deleteat!(P,delete_list)
    return P 
end

#Input data
L = 4
data_label = "ovarian.txt"

#-----main function starts----

transitions = possible_transitions(L)
n_partners = transitions.n_partners
cumulative_partners = transitions.cumulative_partners
partners = transitions.partners
edges_list = edges(L,n_partners,cumulative_partners,partners)
start = edges_list.start
dest = edges_list.dest
var = define_variables(L, start, dest)
a = var.a
b = var.b
N = read_data(L,data_label)
As = A_polynomials(L,var.R,a,b,N,start,dest,var.alist,var.blist)

P = polynomial_system(L,var.R,a,b,As,start,dest,var.alist,var.blist)

P_red = remove_variables(n_partners,start,dest,L,a,var.alist,b,var.blist,P.P_list,P.Pb_list,P.index_list)
P_final = remove_polynomials(P_red,L)

I = ideal(P_final)
#------- Evaluation for ovarian example HyperLAU

G = groebner_basis(I)

values_HyperLAU = [0.5024,0.1398,0.2157,0.142,0.4973,0.0468,0.4558,0.21,0.052,0.738,0.9999,0,0.2425,0.5212,0.2363,0.2365,0.7635,0.7596,0.2404,0.8707,0,0.1292,0,0.9999,0.0005,0.9995,0.558,0.442]
values_HyperTraPS = [0.4684,0.07,0.2561,0.1984,0.6125,0.081,0.3066,0.0092,0.0006,0.9903,0.9538,0.0462,0.2878,0.703,0.0092,0.309,0.691,0.1093,0.8907,0.8249,0.0733,0.1019,0.908,0.092,0.0115,0.9885,0.9771,0.0229]
subs_dict = QQ.(rationalize.(values_HyperTraPS))
P_1 = evaluate(P_final[1],a[1:28],subs_dict)
P_2 = evaluate(P_final[2],a[1:28],subs_dict)
P_3 = evaluate(P_final[3],a[1:28],subs_dict)
P_4 = evaluate(P_final[4],a[1:28],subs_dict)
P_5 = evaluate(P_final[5],a[1:28],subs_dict)
P_6 = evaluate(P_final[6],a[1:28],subs_dict)
P_7 = evaluate(P_final[7],a[1:28],subs_dict)
P_8 = evaluate(P_final[8],a[1:28],subs_dict)
P_9 = evaluate(P_final[9],a[1:28],subs_dict)
P_10 = evaluate(P_final[10],a[1:28],subs_dict)
P_11 = evaluate(P_final[11],a[1:28],subs_dict)
P_12 = evaluate(P_final[12],a[1:28],subs_dict)
P_13 = evaluate(P_final[13],a[1:28],subs_dict)
P_14 = evaluate(P_final[14],a[1:28],subs_dict)
P_15 = evaluate(P_final[15],a[1:28],subs_dict)
P_16 = evaluate(P_final[16],a[1:28],subs_dict)
P_17 = evaluate(P_final[17],a[1:28],subs_dict)
P_18 = evaluate(P_final[18],a[1:28],subs_dict)
P_19 = evaluate(P_final[19],a[1:28],subs_dict)
P_20 = evaluate(P_final[20],a[1:28],subs_dict)
P_21 = evaluate(P_final[21],a[1:28],subs_dict)
P_22 = evaluate(P_final[22],a[1:28],subs_dict)
P_23 = evaluate(P_final[23],a[1:28],subs_dict)
P_24 = evaluate(P_final[24],a[1:28],subs_dict)
P_25 = evaluate(P_final[25],a[1:28],subs_dict)
P_26 = evaluate(P_final[26],a[1:28],subs_dict)
P_27 = evaluate(P_final[27],a[1:28],subs_dict)
P_28 = evaluate(P_final[28],a[1:28],subs_dict)

b_HyperLAU = [0.5206,0.1785,0.2872,1,0.1418,0.7128,0.5543,0.0058,0.0717,0.8582,0.8215,0.0426,0.5937,0.1945,0.4005,0.0515,0.9283,0,0.4794,0.0025,0.3353,0.0014,0.1104,0.3992,0.8041,0.9549,0.2202,0.3291]
b_HyperTraPS = [0.7154,1,0.193,0.1235,0.1005,0.807,0.8829,0.2675,0.0484,0.8995,0,0.0219,0.2037,0.0069,0.5288,0.0868,0.9516,0.8765,0.2846,0.9163,0.0579,0.0001,0.0592,0.3293,0.993,0.0618,0.1667,0.4172]
subs_dict_2 = QQ.(rationalize.(b_HyperTraPS))

P1 = evaluate(P_1,b[1:28],subs_dict_2)
P2 = evaluate(P_2,b[1:28],subs_dict_2)
P3 = evaluate(P_3,b[1:28],subs_dict_2)
P4 = evaluate(P_4,b[1:28],subs_dict_2)
P5 = evaluate(P_5,b[1:28],subs_dict_2)
P6 = evaluate(P_6,b[1:28],subs_dict_2)
P7 = evaluate(P_7,b[1:28],subs_dict_2)
P8 = evaluate(P_8,b[1:28],subs_dict_2)
P9 = evaluate(P_9,b[1:28],subs_dict_2)
P10 = evaluate(P_10,b[1:28],subs_dict_2)
P11 = evaluate(P_11,b[1:28],subs_dict_2)
P13 = evaluate(P_13,b[1:28],subs_dict_2)
P14 = evaluate(P_14,b[1:28],subs_dict_2)
P15 = evaluate(P_15,b[1:28],subs_dict_2)
P16 = evaluate(P_16,b[1:28],subs_dict_2)
P17 = evaluate(P_17,b[1:28],subs_dict_2)
P18 = evaluate(P_18,b[1:28],subs_dict_2)
P19 = evaluate(P_19,b[1:28],subs_dict_2)
P20 = evaluate(P_20,b[1:28],subs_dict_2)
P21 = evaluate(P_21,b[1:28],subs_dict_2)
P22 = evaluate(P_22,b[1:28],subs_dict_2)
P23 = evaluate(P_23,b[1:28],subs_dict_2)
P24 = evaluate(P_24,b[1:28],subs_dict_2)
P25 = evaluate(P_25,b[1:28],subs_dict_2)
P26 = evaluate(P_26,b[1:28],subs_dict_2)
P27 = evaluate(P_27,b[1:28],subs_dict_2)
P28 = evaluate(P_28,b[1:28],subs_dict_2)
