# HyperALG
An Algebraic Approach to Evolutionary Accumulation Models: Creating a set of polynomial equations for given cross-sectional data, whose solution set describes all possible transition parameters to that system. 

<img width="740" height="514" alt="overview" src="https://github.com/user-attachments/assets/1b6d4a11-48e7-4f98-a71f-787f54dd54dc" />

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
