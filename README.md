# Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range

Code repository for Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range.

The following is the list of priority assignment algorithms and schedulability tests implemented in this repository:

## Priority Assignment Methods
- Deadline Monotonic Priority Ordering (DMPO)
- Deadline Minus Computation Monotonic (D-CMPO)
- Deadline Minus kC (DkC): k is a function of the number of processors in the system
- Audsley's Optimal Priority Assignment (OPA): implementation of N.C. Audsley, _On priority assignment in fixed priority scheduling_, https://doi.org/10.1016/S0020-0190(00)00165-4.
- MITER-based Optimization Method: implementation of Y. Zhao and H. Zeng, _The Concept of Response Time Estimation Range for Optimizing Systems Scheduled with Fixed Priority_, https://doi.org/10.1109/RTAS.2018.00036. This method is referred to as MUTER in the code, as per the name given in the RTSS 2018 paper. 
- Hybrid Priority Assignment using MITER (HP-MITER): Hybrid Priority Assignment method using a combination of the MITER-based method and DkC.


## Schedulability Analyses
- Deadline Analysis with Limited Carry-in (DA-LC): implementation of the method defined in Davis, R.I., Burns, A., _Improved priority assignment for global fixed priority pre-emptive scheduling in multiprocessor real-time systems_, https://doi.org/10.1007/s11241-010-9106-5.
- GSYY: implementation of the test defined in N. Guan, M. Stigge, W. Yi and G. Yu, _New Response Time Bounds for Fixed Priority Multiprocessor Scheduling_, https://doi.org/10.1109/RTSS.2009.11.
- ZLL: implementation of the test defined in Q. Zhou, G. Li and J. Li, _Improved Carry-in Workload Estimation for Global Multiprocessor Scheduling_, https://doi.org/10.1109/TPDS.2017.2679195.


## **IMPORTANT NOTE:** 
- The MITER algorithm uses the CPLEX optimization solver, however as it is proprietary, it has not been included in this repository.
- HP-MITER requires the implementation of the EPE-ZLL test (Quan Zhou, Jianjun Li, and Guohui Li, _Excluding Parallel Execution to Improve Global Fixed Priority Response Time Analysis_, https://doi.org/10.1145/3477035). We obtained this implementation from the authors. However, as we do not yet have permission to post their code publicly, this implementation is also not included in this repository. 


## How to run *main.cpp*
- Uncomment one of the existing *Experiment()* calls in the *main()* function in ***main.cpp*** or write a new one. The arguments of the *Experiment()* function are:
  - Lower Bound of Utilization
  - Upper Bound of Utilization
  - Number of Processors
  - Number of Tasks in each Taskset
  - Number of points in the graph between the lower and upper bounds
  - `true` for constrained deadline and `false` for implicit deadline
  - Name of the folder to be created with the test tasksets

- Change the filename in the definition of the variable **_axOutput_**. This is the file where the consolidated results will be written. 

- To build the code, run
```
$ make exec
```

- To run the code, the command is
```
$ ./GlobalFixedPriorityScheduling
```

- After the program terminates successfully, the file named in **_axOutput_** will be created in the root directory. 
- This file can be parsed with a .py script to obtain graphs.


## Acknowledgement

- We thank the authors of _Excluding Parallel Execution to Improve Global Fixed Priority Response Time Analysis_, https://doi.org/10.1145/3477035 for providing the implementation of their work.
 
- We acknowledge Advanced Research Computing at Virginia Tech for providing computational resources and technical support that have contributed to the testing of this work. URL: https://arc.vt.edu/.

## References
\[1] 
