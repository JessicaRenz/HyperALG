#!/bin/bash

#OUTPUT_FILE = "optimization_results.csv"

#Uniform strategies
for a_val in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
    for b_val in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
        a=()
        b=()
        for i in {1..28}; do
            a+=($a_val)
            b+=($b_val)
        done
        a_string=$(IFS=,; echo "${a[*]}")
        b_string=$(IFS=,; echo "${b[*]}")
        julia ./optimizing_a_b.jl $a_string $b_string 
    done
done

#Random strategies
for i in 42 1337 2024 9876 5555 123456 777 31415 99999 271828; do
    a=()
    b=()
    for j in {1..28}; do
        a+=($(awk -v seed=$i 'BEGIN{srand(seed); print rand()}'))
        b+=($(awk -v seed=$i 'BEGIN{srand(seed); print rand()}'))
    done
    a_string=$(IFS=,; echo "${a[*]}")
    b_string=$(IFS=,; echo "${b[*]}")
    julia ./optimizing_a_b.jl $a_string $b_string 
done

# Ranges
for start in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
    for end in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
        a=()
        b=()
        b_alt=()
        for i in {0..27}; do
            val_a=$(awk -v start=$start -v end=$end -v i=$i -v n=27 'BEGIN{print start + i*(end - start)/n}')
            a+=($val_a)
            val_b_alt=$(awk -v start=$end -v end=$start -v i=$i -v n=27 'BEGIN{print start +i*(end - start)/n}')
            b_alt+=($val_b_alt)
        done
        a_string=$(IFS=,; echo "${a[*]}")
        b_string=$(IFS=,; echo "${b[*]}")
        b_alt_string=$(IFS=,; echo "${b_alt[*]}")
        julia ./optimizing_a_b.jl $a_string $a_string 
        julia ./optimizing_a_b.jl $a_string $b_alt_string 
    done
done
