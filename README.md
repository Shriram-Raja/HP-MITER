# Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range

Code repository for Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range.

The following is the list of priority assignment algorithms and schedulability tests implemented in this repository:

### Priority Assignment Methods
- Deadline Monotonic Priority Ordering (DMPO)
- Deadline Minus Computation Monotonic (D-CMPO)
- Deadline Minus kC (DkC): k is a function of the number of processors in the system
- Audsley's Optimal Priority Assignment (OPA)
- MITER-based Optimization Method (MITER): Implementation of _**The Concept of Response Time Estimation Range for Optimizing Systems Scheduled with Fixed Priority**_, RTSS 2018
- Hybrid Priority Assignment using MITER (HP-MITER): Hybrid Priority Assignment method using a combination of MITER and DkC

### Schedulability Analyses
- Deadline Analysis with Limited Carry-in (DA-LC): implementation of Davis, R.I., Burns, A. _Improved priority assignment for global fixed priority pre-emptive scheduling in multiprocessor real-time systems_ Real-Time Syst 47, 1–40 (2011). https://doi.org/10.1007/s11241-010-9106-5
- GSYY: implementation of N. Guan, M. Stigge, W. Yi and G. Yu, "_New Response Time Bounds for Fixed Priority Multiprocessor Scheduling_", 2009 30th IEEE Real-Time Systems Symposium, Washington, DC, USA, 2009, pp. 387-397, doi: 10.1109/RTSS.2009.11.
- ZLL: Implementation of Q. Zhou, G. Li and J. Li, "_Improved Carry-in Workload Estimation for Global Multiprocessor Scheduling_", in IEEE Transactions on Parallel and Distributed Systems, vol. 28, no. 9, pp. 2527-2538, 1 Sept. 2017, doi: 10.1109/TPDS.2017.2679195.

**IMPORTANT NOTE:** 
- The MITER algorithm uses the CPLEX optimization solver, however as it is proprietary, it has not been included in this repository.
- HP-MITER requires the implementation of EPE-ZLL (Quan Zhou, Jianjun Li, and Guohui Li. 2021. Excluding Parallel Execution to Improve Global Fixed Priority Response Time Analysis. ACM Trans. Embed. Comput. Syst. 20, 5s, Article 104 (October 2021), 24 pages. https://doi.org/10.1145/3477035). We obtained this implementation from the authors. As we have not obtained permission to make their code public, this is also not included in this repository. 

## How to run *main.cpp*

- Uncomment one of the existing *Experiment()* calls in the *main()* function in ***main.cpp*** or write a new one. The arguments of the *Experiment()* function are:
  - Lower Bound of Utilization
  - Upper Bound of Utilization
  - Number of Processors
  - Number of Tasks in each Taskset
  - Number of points in the graph between the lower and upper bounds
  - `true` for constrained deadline and `false` for implicit deadline
  - Name of the folder to be created with the test tasksets

- Change the filename in Line 49 of the file. This is the file where the consolidated results will be written into. 

- To build the code, run
```
$ make exec
```

- To run the code, the command is
```
$ ./GlobalFixedPriorityScheduling
```

- After the run is complete, the file named in Line 49 will be created in the root directory. 
- This file can be parsed with a .py script to obtain graphs. 

## Modifying Priority Orderings (POs)

- Open ***RandomSystemsTest.cpp***.
- *TestUSweep()*:
  - The function *TestUSweep()*, specifically the part from Line 101 to Line 203, is modified to change which POs are computed.
  - The part of the function from Line 204 to the end of the file should be modified to change which POs are written into the Result.txt file for each taskset.
- *MarshallUSweepResult()*:
  - The function *MarshallUSweepResult()* in the same file is reads all the Result.txt files and compiles them into a single file. 
  - If the number of POs used has been modified in the *TestUSweep()* function, the number in the condition of the inner `for` loops should be modified.
