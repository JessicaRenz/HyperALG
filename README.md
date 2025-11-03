# HyperALG
An Algebraic Approach to Evolutionary Accumulation Models: Creating a set of polynomial equations for given cross-sectional data, whose solution set describes all possible transition parameters to that system. 

<img width="774" height="543" alt="overview" src="https://github.com/user-attachments/assets/83fc42f2-4593-4756-b530-1432266eaf0f" />


## Requirements
For running HyperALG, you need the ability to run Julia code, and having the packages `Oscar.jl` and `DelimitedFiles.jl`.

## Running HyperALG
HyperALG can be runned directly from the command line by specifying the data set as an input parameter:
```
julia HyperALG.jl [dataset.txt]
```
 HyperLAU expects the data as a text file, the ending has to be included in the input parameter. Make sure to save the dataset in the same directory as the `HyperALG.jl` file.

 Some functions in this code originally come from HyperHMM https://github.com/StochasticBiology/hypercube-hmm/tree/main and were transferred to Julia. The functions to which this applies are marked in the source code.

 The input is expected to be a list of binary strings (cross-sectional data), for example:
 ```
001
001
100
101
```

 ## Output
HyperALG outputs the three text files `polynomials_[dataset.txt]`, `a_variables_[dataset.txt]` and `b_variables_[dataset.txt]`.
- **polynomials_[...].txt** In this output you can find the final set of polynomials describing the evolutionary process given the data. Setting all of them `= 0` gives you the set of polynomial equations that the transition parameters need to fulfill. The variables `a[i]` describe the forward transitions and the variables `b[i]` the backwards proportions.
- **a_variables_[...].txt** This file contains a table that shows the assignment of the `a[i]` to the edges in the hypercube. The first column contains the name of the variable, the second and third contain the start node and the end node of the corresponding edge, respectively.
- **b_variables_[...].txt** This file contains a table that shows the assignment of the `b[i]` to the edges in the hypercube. The first column contains the name of the variable, the second and third contain the start node and the end node of the corresponding edge, respectively. 

## Case Studies
In this folder you can find the data and code we used for our case studies in the corresponding paper. 

### toy.txt
This file contains the data to the small toy example for `L=3` in the paper.

### ovarian.txt
This file contains the data we used for the cancer case study for `L=4` in the paper. The features describe the presence or absence of chromosomal abberations in 87 samples of ovarian cancer patients: Feature 1: `8q+`, Feature 2: `3q+`, Feature 3: `5q-` and Feature 4: `4q-`. These labels indicate the index of the chromosome where the mutation occurs, as well as the chromosomal arm (p or q). A + indicates addition while - denotes a loss.

### optimizing_a_b.jl
For running this Julia code, the package `Optim.jl` needs to be installed. It takes the 28 polynomials we get when running `HyperALG.jl` on the dataset `ovarian.txt` and creates the sum of the squares of these. Minimizing the obtained function delivers us an estimated value for all veriables `a[i]` and `b[i]`. 

### Ovarian_L4.jl
For running this Julia code, the package `Oscar.jl` needs to be installed. It contains the HyperALG code itself as well as the evaluation of the polynomials obtained at the values provided by HyperLAU (https://github.com/JessicaRenz/HyperLAU) and HyperTraPS (https://github.com/StochasticBiology/hypertraps-ct) for the `a[i]` variables. The results are also based on the dataset `ovarian.txt`. In a second step, also the values for the `b[i]` variables are substituted. These are obtained by an minimization procedure exactly as in `optimizing_a_b.jl`. 

## References
Data source:

Knutsen, T., Gobu, V., Knaus, R., Padilla‐Nash, H., Augustus, M., Strausberg, R.L., Kirsch, I.R., Sirotkin, K. and Ried, T., 2005. The interactive online SKY/M‐FISH & CGH database and the Entrez cancer chromosomes search database: linkage of chromosomal aberrations with the genome sequence. Genes, Chromosomes and Cancer, 44(1), pp.52-64.
