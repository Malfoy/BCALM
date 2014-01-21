bcalm
=====

de Bruijn CompAction in Low Memory

Possible use :

./bcalm input.dot
detect k and do the compaction with m=10 in compacted.dot

./bcalm input.dot output.dot
detect k and do the compaction with m=10 in output.dot

./bcalm input.dot output.dot 8
detect k and do the compaction with m=8 in output.dot



Nota Bene :   
Higher m mean lower memory but the algorithm will NOT work with m>10.   
Also note that to use m>=8 you have to allow 260 open files (`ulimit -n 260`)
and 1110 for m=10 (`ulimit -n 1100`).

Inputs
======

Convert DSK output to bcalm as follows: use the `parse_results` program
provided with DSK with the `--dot` flag, e.g.: 

`./parse_results dsk_output.solid_kmers_binary --dot > input_bcalm.dot`
