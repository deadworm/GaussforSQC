There are 3 input files need to be put in the same directory with "GaussforSQC.exe".
And if there is no wrong, there are 3 output files. I will introduce them one by one.

input files:
1. infile.txt
2. r_gauss.txt
3. set_basis.txt

output files:
4. Fbase.txt
5. Hv.txt
6. Ht.txt

******
1. infile.txt:

1E-05   ->  this is the minimum of the output coefficients

2   ->  this is number of threads used in your calculation

3   ->  this is the number of "atoms" in your calculation

O 6.45  0.00000000  0.00000000  0.00000000    ->  this is the name, nuclear charge, and x-y-z-position of your "atom"
H 1 1.8 0.00000000 0.00000000                     (name can only marked by one letter and should be the same with one in your set_basis.txt file)
H 1  -0.45 1.75 0.00000000  

******
2. r_gauss.txt:

10    ->  this is the number of guass functions used to fit Coulomb interaction
          (for high precision you can write your own file in the same form)
0.02283   ->  this is the constant term

0.99474
1.64996
0.77793
0.84743
1.52939
1.04353
1.39758
2.30849
1.33282
1.39149

1.501
12.74385
0.63876
0.9778
5.78379
0.41838
2.76619
32.48339
0.22675
0.11477

******
3. set_basis.txt:

2   ->  this is the number of "elements" in your calculation

O 2    -> this is the name and number of orbits for your atom

  S 3   ->  this is the name of orbits and number of Gauss functions of this orbits
            (only S,P,D,F can be used for now and those coefficients below can be copied in any basis sets library)
  5.0331513             1.1695961              0.3803890              
  -0.09996723             0.39951283             0.70011547             

  P 3
  5.0331513             1.1695961              0.3803890
  0.15591627       0.60768372       0.39195739       

H 1

  S 3
  3.42525091             0.62391373             0.16885540             
  0.15432897       0.53532814       0.44463454       
  
******
4. Fbase.txt  (this is a file give you the orbits after normalization 
               and can only be useful in analysing componets of those orbits 
               or calculating distribution of electrons)
               
******
5. Hv.txt

nuclear distance
1.8,	2.85044,	  ->  this is the distance between the first and second, second and third "atoms"

nuclear potential energy
7.50374,    ->  this is the repulsive energy of all "atoms"

potential energy matrix
-8.36583,	-0.143254,	-0.183311,	0,	1.07174,	0.928566,	    ->  this is the matrix for SQC one electron potential terms
-0.143254,	-7.24968,	0.0356165,	0,	0.899093,	-0.341199,	
-0.183311,	0.0356165,	-7.2303,	0,	0.00982184,	0.87914,	
0,	0,	0,	-7.09179,	0,	0,	
1.07174,	0.899093,	0.00982184,	0,	-5.17327,	-0.0213473,	
0.928566,	-0.341199,	0.87914,	0,	-0.0213473,	-5.15713,	

kinetic energy matrix
0.808128,	0,	0,	0,	-0.339912,	-0.294615,	  ->  this is the matrix for SQC one electron kinetic energy terms
0,	1.80707,	0,	0,	-0.653684,	0.255078,	
0,	0,	1.80707,	0,	0,	-0.638094,	
0,	0,	0,	1.80707,	0,	0,	
-0.339912,	-0.653684,	0,	0,	1.58085,	-0.130072,	
-0.294615,	0.255078,	-0.638094,	0,	-0.130072,	1.61072,

******
6. Ht.txt   (this file has almost the same form with Hv.txt and keeps SQC two electrons interaction energy terms)




So, with those coefficients for normal form of second quantization Hamiltonian in output files Hv.txt and Ht.txt,
The numerical result of ground state of your system can be get by exact diagonalization or other methods like DMRG.











