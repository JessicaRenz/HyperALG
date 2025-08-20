using Oscar

R, (a04, a02, a15, a26, a46, b15, b13, b26, b37, b57) = polynomial_ring(QQ,["a04","a02","a15","a26","a46","b15","b13","b26","b37","b57"])

D = [1,2,0,3,1,2,0,1]

D0 = D[1]
D1 = D[2]
D2 = D[3]
D3 = D[4]
D4 = D[5]
D5 = D[6]
D6 = D[7]
D7 = D[8]

n = D0 + D1 + D2 + D3 + D4 + D5 + D6 + D7

N0 = rationalize(D0/n)
N1 = rationalize(D1/n)
N2 = rationalize(D2/n)
N3 = rationalize(D3/n)
N4 = rationalize(D4/n)
N5 = rationalize(D5/n)
N6 = rationalize(D6/n)
N7 = rationalize(D7/n)

A4 = N0*a04 + N4 + N5*(1-b15) + N6*(1-b26) + N7*(1-b15)*b57 + N7*(1-b26)*(1-b37-b57)
A2 = N0*a02 + N2 + N3*(1-b13) + N6*b26 + N7*(1-b13)*b37 + N7*b26*(1-b37-b57)
A1 = N0*(1-a04-a02) + N1 + N3*b13 + N5*b15 + N7*b13*b37 + N7*b15*b57
A6 = N0*a02*a26 + N0*a04*a46 + N2*a26 + N4*a46 + N6 + N7*(1-b37-b57)
A5 = N0*(1-a04-a02)*a15 + N0*a04*(1-a46) + N1*a15 + N4*(1-a46) + N5 + N7*b57
A3 = N0*(1-a02-a04)*(1-a15) + N0*a02*(1-a26) + N1*(1-a15) + N2*(1-a26) + N3 + N7*b37

p1 = 1-a04-a02-A1
p2 = a02 - A2
p4 = A1*(1-a15) + A2*(1-a26) - A3
p5 = A1*a15 + A4*(1-a46) - A5
p7 = b15*(A1*a15 + A4*(1-a46)) - A1*a15
p8 = b13*(A1*(1-a15) + A2*(1-a26)) - A1*(1-a15)
p9 = b26*(A2*a26 + A4*a46) - A2*a26
p10 = b37*(A3 + A5 + A6) - A3
p11 = b57*(A3 + A5 + A6) - A5