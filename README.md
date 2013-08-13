bcalm
=====

de Bruijn CompAction in Low Memory

Possible use :

./bcalm   
Launch a small test on a random genome


./bcalm n   
Launch a test on a random genome of n kmers


./bcalm n k   
Same as before but with a specified value for k


./bcalm n k m   
Same as before but with a specified value of the minimiser size (<=10)


./bcalm input output k m    
Real utilisation : read input file and write the nodes in the output file

Nota Bene :   
Higher m mean lower memory but the algorithm will NOT work with m>10.   
Also note that to use m>=8 you have to allow 260 open files (with ulimit -n 260)
and 1024 to use m=10 (with ulimit -n 1024).
