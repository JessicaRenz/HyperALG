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


#------- Evaluation for ovarian example HyperLAU

values_HyperLAU = [0.502,0.14,0.216,0.142,0.497,0.047,0.456,0.21,0.052,0.738,1,0,0.243,0.521,0.236,0.236,0.764,0.76,0.24,0.871,0,0.129,0,1,0.001,0.999,0.558,0.442]
subs_dict = QQ.(rationalize.(values_HyperLAU))
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

# P_26 => b[16] = 43367147//837500000 = 0.052
# P_27 => b[24] = 66748591//167500000 = 0.398
# P_28 => b[27] = 368999127//1675000000 = 0.22

vars_to_subs_1 = [b[16],b[24],b[27]]
vals_to_subs_1 = QQ.([43367147//837500000, 66748591//167500000,368999127//1675000000])
P_1_1 = evaluate(P_1,vars_to_subs_1,vals_to_subs_1)
P_2_1 = evaluate(P_2,vars_to_subs_1,vals_to_subs_1)
P_3_1 = evaluate(P_3,vars_to_subs_1,vals_to_subs_1)
P_4_1 = evaluate(P_4,vars_to_subs_1,vals_to_subs_1)
P_5_1 = evaluate(P_5,vars_to_subs_1,vals_to_subs_1)
P_6_1 = evaluate(P_6,vars_to_subs_1,vals_to_subs_1)
P_7_1 = evaluate(P_7,vars_to_subs_1,vals_to_subs_1)
P_8_1 = evaluate(P_8,vars_to_subs_1,vals_to_subs_1)
P_9_1 = evaluate(P_9,vars_to_subs_1,vals_to_subs_1)
P_10_1 = evaluate(P_10,vars_to_subs_1,vals_to_subs_1)
P_11_1 = evaluate(P_11,vars_to_subs_1,vals_to_subs_1)
P_12_1 = evaluate(P_12,vars_to_subs_1,vals_to_subs_1)
P_13_1 = evaluate(P_13,vars_to_subs_1,vals_to_subs_1)
P_14_1 = evaluate(P_14,vars_to_subs_1,vals_to_subs_1)
P_15_1 = evaluate(P_15,vars_to_subs_1,vals_to_subs_1)
P_16_1 = evaluate(P_16,vars_to_subs_1,vals_to_subs_1)
P_17_1 = evaluate(P_17,vars_to_subs_1,vals_to_subs_1)
P_18_1 = evaluate(P_18,vars_to_subs_1,vals_to_subs_1)
P_19_1 = evaluate(P_19,vars_to_subs_1,vals_to_subs_1)
P_20_1 = evaluate(P_20,vars_to_subs_1,vals_to_subs_1)
P_21_1 = evaluate(P_21,vars_to_subs_1,vals_to_subs_1)
P_22_1 = evaluate(P_22,vars_to_subs_1,vals_to_subs_1)
P_23_1 = evaluate(P_23,vars_to_subs_1,vals_to_subs_1)
P_24_1 = evaluate(P_24,vars_to_subs_1,vals_to_subs_1)
P_25_1 = evaluate(P_25,vars_to_subs_1,vals_to_subs_1)

#P_15 and P_7 => b[8] = 0

vars_to_subs_2 = [b[8]]
vals_to_subs_2 = QQ.([0])
P_1_2 = evaluate(P_1_1,vars_to_subs_2,vals_to_subs_2)
P_2_2 = evaluate(P_2_1,vars_to_subs_2,vals_to_subs_2)
P_3_2 = evaluate(P_3_1,vars_to_subs_2,vals_to_subs_2)
P_4_2 = evaluate(P_4_1,vars_to_subs_2,vals_to_subs_2)
P_5_2 = evaluate(P_5_1,vars_to_subs_2,vals_to_subs_2)
P_6_2 = evaluate(P_6_1,vars_to_subs_2,vals_to_subs_2)
P_7_2 = evaluate(P_7_1,vars_to_subs_2,vals_to_subs_2)
P_8_2 = evaluate(P_8_1,vars_to_subs_2,vals_to_subs_2)
P_9_2 = evaluate(P_9_1,vars_to_subs_2,vals_to_subs_2)
P_10_2 = evaluate(P_10_1,vars_to_subs_2,vals_to_subs_2)
P_11_2 = evaluate(P_11_1,vars_to_subs_2,vals_to_subs_2)
P_12_2 = evaluate(P_12_1,vars_to_subs_2,vals_to_subs_2)
P_13_2 = evaluate(P_13_1,vars_to_subs_2,vals_to_subs_2)
P_14_2 = evaluate(P_14_1,vars_to_subs_2,vals_to_subs_2)
P_15_2 = evaluate(P_15_1,vars_to_subs_2,vals_to_subs_2)
P_16_2 = evaluate(P_16_1,vars_to_subs_2,vals_to_subs_2)
P_17_2 = evaluate(P_17_1,vars_to_subs_2,vals_to_subs_2)
P_18_2 = evaluate(P_18_1,vars_to_subs_2,vals_to_subs_2)
P_19_2 = evaluate(P_19_1,vars_to_subs_2,vals_to_subs_2)
P_20_2 = evaluate(P_20_1,vars_to_subs_2,vals_to_subs_2)
P_21_2 = evaluate(P_21_1,vars_to_subs_2,vals_to_subs_2)
P_22_2 = evaluate(P_22_1,vars_to_subs_2,vals_to_subs_2)
P_23_2 = evaluate(P_23_1,vars_to_subs_2,vals_to_subs_2)
P_24_2 = evaluate(P_24_1,vars_to_subs_2,vals_to_subs_2)
P_25_2 = evaluate(P_25_1,vars_to_subs_2,vals_to_subs_2)

# P_10 and P_19 => b[7] = 25919000//46748591 = 0.554

vars_to_subs_3 = [b[7]]
vals_to_subs_3 = QQ.([25919000//46748591])
P_1_3 = evaluate(P_1_2,vars_to_subs_3,vals_to_subs_3)
P_2_3 = evaluate(P_2_2,vars_to_subs_3,vals_to_subs_3)
P_3_3 = evaluate(P_3_2,vars_to_subs_3,vals_to_subs_3)
P_4_3 = evaluate(P_4_2,vars_to_subs_3,vals_to_subs_3)
P_5_3 = evaluate(P_5_2,vars_to_subs_3,vals_to_subs_3)
P_6_3 = evaluate(P_6_2,vars_to_subs_3,vals_to_subs_3)
P_7_3 = evaluate(P_7_2,vars_to_subs_3,vals_to_subs_3)
P_8_3 = evaluate(P_8_2,vars_to_subs_3,vals_to_subs_3)
P_9_3 = evaluate(P_9_2,vars_to_subs_3,vals_to_subs_3)
P_10_3 = evaluate(P_10_2,vars_to_subs_3,vals_to_subs_3)
P_11_3 = evaluate(P_11_2,vars_to_subs_3,vals_to_subs_3)
P_12_3 = evaluate(P_12_2,vars_to_subs_3,vals_to_subs_3)
P_13_3 = evaluate(P_13_2,vars_to_subs_3,vals_to_subs_3)
P_14_3 = evaluate(P_14_2,vars_to_subs_3,vals_to_subs_3)
P_16_3 = evaluate(P_16_2,vars_to_subs_3,vals_to_subs_3)
P_17_3 = evaluate(P_17_2,vars_to_subs_3,vals_to_subs_3)
P_18_3 = evaluate(P_18_2,vars_to_subs_3,vals_to_subs_3)
P_19_3 = evaluate(P_19_2,vars_to_subs_3,vals_to_subs_3)
P_20_3 = evaluate(P_20_2,vars_to_subs_3,vals_to_subs_3)
P_21_3 = evaluate(P_21_2,vars_to_subs_3,vals_to_subs_3)
P_22_3 = evaluate(P_22_2,vars_to_subs_3,vals_to_subs_3)
P_23_3 = evaluate(P_23_2,vars_to_subs_3,vals_to_subs_3)
P_24_3 = evaluate(P_24_2,vars_to_subs_3,vals_to_subs_3)
P_25_3 = evaluate(P_25_2,vars_to_subs_3,vals_to_subs_3)

# b[12] = 502301107673322//11455033051430669 = 0.044
# b[13] = 6917990836007349//11455033051430669 = 0.604
# b[14] =  2483462154057092052168420//12076819908879437120437561 = 0.206
# b[20] = 0
# b[21] = 124088503869084686272548421616073450//363723236920488507315506410982483269 = 0.341

vars_to_subs_4 = [b[12],b[13],b[14],b[20],b[21]]
vals_to_subs_4 = QQ.(rationalize.([0.044,0.604,0.206,0,0.341]))
P_1_4 = evaluate(P_1_3, vars_to_subs_4,vals_to_subs_4)
P_2_4 = evaluate(P_2_3, vars_to_subs_4,vals_to_subs_4)
P_3_4 = evaluate(P_3_3, vars_to_subs_4,vals_to_subs_4)
P_4_4 = evaluate(P_4_3, vars_to_subs_4,vals_to_subs_4)
P_5_4 = evaluate(P_5_3, vars_to_subs_4,vals_to_subs_4)
P_6_4 = evaluate(P_6_3, vars_to_subs_4,vals_to_subs_4)
P_7_4 = evaluate(P_7_3, vars_to_subs_4,vals_to_subs_4)
P_8_4 = evaluate(P_8_3, vars_to_subs_4,vals_to_subs_4)
P_9_4 = evaluate(P_9_3, vars_to_subs_4,vals_to_subs_4)
P_10_4 = evaluate(P_10_3, vars_to_subs_4,vals_to_subs_4)
P_11_4 = evaluate(P_11_3, vars_to_subs_4,vals_to_subs_4)
P_12_4 = evaluate(P_12_3, vars_to_subs_4,vals_to_subs_4)
P_13_4 = evaluate(P_13_3, vars_to_subs_4,vals_to_subs_4)
P_14_4 = evaluate(P_14_3, vars_to_subs_4,vals_to_subs_4)
P_16_4 = evaluate(P_16_3, vars_to_subs_4,vals_to_subs_4)
P_17_4 = evaluate(P_17_3, vars_to_subs_4,vals_to_subs_4)
P_18_4 = evaluate(P_18_3, vars_to_subs_4,vals_to_subs_4)
P_19_4 = evaluate(P_19_3, vars_to_subs_4,vals_to_subs_4)
P_20_4 = evaluate(P_20_3, vars_to_subs_4,vals_to_subs_4)
P_21_4 = evaluate(P_21_3, vars_to_subs_4,vals_to_subs_4)
P_22_4 = evaluate(P_22_3, vars_to_subs_4,vals_to_subs_4)
P_23_4 = evaluate(P_23_3, vars_to_subs_4,vals_to_subs_4)
P_24_4 = evaluate(P_24_3, vars_to_subs_4,vals_to_subs_4)
P_25_4 = evaluate(P_25_3, vars_to_subs_4,vals_to_subs_4)

# b[22] = 0.0001

vars_to_subs_5 = [b[22]]
vals_to_subs_5 = QQ.(rationalize.([0.0001]))
P_1_5 = evaluate(P_1_4,vars_to_subs_5,vals_to_subs_5)
P_2_5 = evaluate(P_2_4,vars_to_subs_5,vals_to_subs_5)
P_3_5 = evaluate(P_3_4,vars_to_subs_5,vals_to_subs_5)
P_4_5 = evaluate(P_4_4,vars_to_subs_5,vals_to_subs_5)
P_5_5 = evaluate(P_5_4,vars_to_subs_5,vals_to_subs_5)
P_6_5 = evaluate(P_6_4,vars_to_subs_5,vals_to_subs_5)
P_8_5 = evaluate(P_8_4,vars_to_subs_5,vals_to_subs_5)
P_9_5 = evaluate(P_9_4,vars_to_subs_5,vals_to_subs_5)
P_10_5 = evaluate(P_10_4,vars_to_subs_5,vals_to_subs_5)
P_11_5 = evaluate(P_11_4,vars_to_subs_5,vals_to_subs_5)
P_12_5 = evaluate(P_12_4,vars_to_subs_5,vals_to_subs_5)
P_13_5 = evaluate(P_13_4,vars_to_subs_5,vals_to_subs_5)
P_14_5 = evaluate(P_14_4,vars_to_subs_5,vals_to_subs_5)
P_17_5 = evaluate(P_17_4,vars_to_subs_5,vals_to_subs_5)
P_18_5 = evaluate(P_18_4,vars_to_subs_5,vals_to_subs_5)
P_19_5 = evaluate(P_19_4,vars_to_subs_5,vals_to_subs_5)
P_20_5 = evaluate(P_20_4,vars_to_subs_5,vals_to_subs_5)
P_21_5 = evaluate(P_21_4,vars_to_subs_5,vals_to_subs_5)
P_22_5 = evaluate(P_22_4,vars_to_subs_5,vals_to_subs_5)
P_24_5 = evaluate(P_24_4,vars_to_subs_5,vals_to_subs_5)
P_25_5 = evaluate(P_25_4,vars_to_subs_5,vals_to_subs_5)

#Test b[5]=1,b[4]=0
vars_to_subs_test = [b[4],b[3],b[5]]
vals_to_subs_test = QQ.(rationalize.([0,0.406,1]))
P_1_test = evaluate(P_1_5,vars_to_subs_test,vals_to_subs_test)
P_2_test = evaluate(P_2_5,vars_to_subs_test,vals_to_subs_test)
P_3_test = evaluate(P_3_5,vars_to_subs_test,vals_to_subs_test)
P_4_test = evaluate(P_4_5,vars_to_subs_test,vals_to_subs_test)
P_5_test = evaluate(P_5_5,vars_to_subs_test,vals_to_subs_test)
P_6_test = evaluate(P_6_5,vars_to_subs_test,vals_to_subs_test)
P_8_test = evaluate(P_8_5,vars_to_subs_test,vals_to_subs_test)
P_9_test = evaluate(P_9_5,vars_to_subs_test,vals_to_subs_test)
P_12_test = evaluate(P_12_5,vars_to_subs_test,vals_to_subs_test)
P_13_test = evaluate(P_13_5,vars_to_subs_test,vals_to_subs_test)
P_14_test = evaluate(P_14_5,vars_to_subs_test,vals_to_subs_test)
P_17_test = evaluate(P_17_5,vars_to_subs_test,vals_to_subs_test)
P_18_test = evaluate(P_18_5,vars_to_subs_test,vals_to_subs_test)
P_21_test = evaluate(P_21_5,vars_to_subs_test,vals_to_subs_test)

vars_bs = [b[1],b[2],b[3],b[16],b[24],b[27],b[8],b[7],b[12],b[13],b[14],b[20],b[21],b[22],b[5],b[4],b[9]]
vals_bs = QQ.(rationalize.([0.1,1,0.428,0.052,0.398,0.22,0,0.495,0.044,0.604,0.206,0.001,0.341,0.001,0,1,0.06]))
P1 = evaluate(P_1,vars_bs,vals_bs)
P2 = evaluate(P_2,vars_bs,vals_bs)
P3 = evaluate(P_3,vars_bs,vals_bs)
P4 = evaluate(P_4,vars_bs,vals_bs)
P5 = evaluate(P_5,vars_bs,vals_bs)
P6 = evaluate(P_6,vars_bs,vals_bs)
P7 = evaluate(P_7,vars_bs,vals_bs)
P8 = evaluate(P_8,vars_bs,vals_bs)
P9 = evaluate(P_9,vars_bs,vals_bs)
P10 = evaluate(P_10,vars_bs,vals_bs)
P11 = evaluate(P_11,vars_bs,vals_bs)
P12 = evaluate(P_12,vars_bs,vals_bs)
P13 = evaluate(P_13,vars_bs,vals_bs)
P14 = evaluate(P_14,vars_bs,vals_bs)
P15 = evaluate(P_15,vars_bs,vals_bs)
P16 = evaluate(P_16,vars_bs,vals_bs)
P17 = evaluate(P_17,vars_bs,vals_bs)
P18 = evaluate(P_18,vars_bs,vals_bs)
P19 = evaluate(P_19,vars_bs,vals_bs)
P20 = evaluate(P_20,vars_bs,vals_bs)
P21 = evaluate(P_21,vars_bs,vals_bs)
P22 = evaluate(P_22,vars_bs,vals_bs)
P23 = evaluate(P_23,vars_bs,vals_bs)
P24 = evaluate(P_24,vars_bs,vals_bs)
P25 = evaluate(P_25,vars_bs,vals_bs)
P26 = evaluate(P_26,vars_bs,vals_bs)
P27 = evaluate(P_27,vars_bs,vals_bs)
P28 = evaluate(P_28,vars_bs,vals_bs)
