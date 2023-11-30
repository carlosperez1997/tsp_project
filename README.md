# tsp_project

# To install graph-tool in MAC M1 (this is a hell!)
!conda install -c conda-forge graph-tool

**'pcb442.tsp'** # it is symmetrical
Objective value: 50783.547 (SAME)
Execution Time: 575.70 s (10 min)

This one has avoided: 706 subtours of size 2 and 126 subtours of size 3

**'u574.tsp'** # it is symmetrical
Objective value: 36934.771 (not same) 
Execution Time: 600.05 (10 min)

This one has avoided: 918 subtours of size 2 and 206 subtours of size 3

**'gr431.tsp'** # symmetrical (geographical)
Objective value: 172042.98 (not optimal, I timed out)
Execution Time: 481.06
This is interesting because we reach a solution. 

Avoided: 823 subtors of 2 and 147 subtors of 3

Question 2:

Does adding constraints 2 and 3, help us reach the solution faster?? Or adding unnecessary constraints takes more time?

Question 3:

Does only adding useful constraints, help us reach the solution faster? Are we able to know which constraints are never reached or violated?