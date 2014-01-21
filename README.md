bcalm
=====

de Bruijn CompAction in Low Memory

Possible use :

./bcalm input.dot
compact input.dot with l=10 in compacted.dot

./bcalm input.dot output.dot
compact input.dot with l=10 in output.dot

./bcalm input.dot output.dot 8
compact input.dot with l=8 in output.dot



Nota Bene :   
Higher l mean lower memory but the algorithm will NOT work with l>10.   
Also note that to use l>=8 you have to allow 260 open files (`ulimit -n 260`)
and 1110 for l=10 (`ulimit -n 1100`).

Inputs
======

Convert DSK output to bcalm as follows: use the `parse_results` program
provided with DSK with the `--dot` flag, e.g.: 

`./parse_results dsk_output.solid_kmers_binary --dot > input_bcalm.dot`
