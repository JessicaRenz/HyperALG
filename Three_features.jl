using Symbolics

function polynomials_L3(D)

    @variables w02 w04 w15 w26 w46 b23 b45 b46 b37 b57 D0 D1 D2 D3 D4 D5 D6 D7 

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

    A4 = N0*w04 + N4 + N6*b46 + N5*b45 + N7*b45*b57 + N7*b46*(1-b37-b57)
    A2 = N0*w02 + N2 + N6 - N6*b46 + N3*b23+ N7*b23*b37 + N7*b46*(1-b37-b57)
    A1 = N0 - N0*w02 - N0*w04 + N1 + N5 - N5*b45 + N3 - N3*b23 + N7*(1-b23)*b37 + N7*(1-b45)*b57
    A6 = N0*w04*w46 + N0*w02*w26 + N4*w46 + N2*w26 + N6 + N7*(1-b37-b57)
    A5 = N0*w04 - N0*w04*w46 + N0*w15 - N0*w02*w15 - N0*w04*w15 + N4 - N4*w46 + N1*w15 + N5 + N7*b57
    A3 = N0 - N0*w02*w26 - N0*w15 - N0*w04 + N0*w04*w15 + N0*w02*w15 + N2 - N2*w26 + N1 - N1*w15 + N3 + N7*b37


    p1 = w04 - A4
    p2 = w02 - A2
    p3 = (A4*w46 + A2*w26) - A6
    p4 = (A4 - A4*w46 + A1*w15) - A5
    p5 = b46*(A4*w46 + A2*w26) - A4*w46
    p6 = b23*(A2*(1-w26) + A1*(1-w15)) - A2*(1-w26)
    p7 = b45*(A1*w15 + A4*(1-w46)) - A4*(1-w46)
    p8 = (A3 + A5 + A6)*b37 - A3
    p9 = (A3 + A5 + A6)*b57 - A5

    equ = [p1,p2,p3,p4,p5,p6,p7,p8,p9]
    return equ 
end 


D = [0,2,0,0,1,5,0,0]
polynomials_L3(D)
