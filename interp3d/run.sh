rm *.sch log history2.txt
mpirun -np 8 nek5000 >> log
