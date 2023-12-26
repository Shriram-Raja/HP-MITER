# Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range

Code repository for Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range.

The following is the list of priority assignment algorithms and schedulability tests implemented in this repository:

## Priority Assignment Methods
- Deadline Monotonic Priority Ordering (DMPO)
- Deadline Minus Computation Monotonic (D-CMPO)
- Deadline Minus kC (DkC): k is a function of the number of processors in the system
- Audsley's Optimal Priority Assignment (OPA): implementation of [[1]](#References).
- MITER-based Optimization Method: implementation of [[2]](#References). This method is called MUTER in the code, as per the name given in the RTSS 2018 paper. 
- Hybrid Priority Assignment using MITER (HP-MITER): Hybrid Priority Assignment method using a combination of the MITER-based method and DkC. This is called HP-MUTER in the code. 

Refer to _MultiProcessorSystem.cpp_

**IMPORTANT NOTE:** The MITER algorithm uses the CPLEX optimization solver, however as it is proprietary, it has not been included in this repository.

## Schedulability Analyses
- Deadline Analysis with Limited Carry-in (DA-LC): implementation of the method defined in [[3]](#References).
- GSYY: implementation of the test defined in [[4]](#References).
- ZLL: implementation of the test defined in [[5]](#References).
- EPE-ZLL: implementation of the test defined in [[6]](#References). This is called GSYY2 in the code.

Refer to _SchedulabilityAnalysis.cpp_


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

- We thank the authors of [[6]](#References) for providing the implementation of their work.
 
- We acknowledge Advanced Research Computing at Virginia Tech for providing computational resources and technical support that have contributed to the testing of this work. URL: https://arc.vt.edu/.

## References

\[1] N.C. Audsley. _On Priority Assignment in Fixed Priority Scheduling_. Information Processing Letters, 79(1):39–44, 2001, https://doi.org/10.1016/S0020-0190(00)00165-4.
<br /><br />
\[2] Y. Zhao and H. Zeng. _The Concept of Response Time Estimation Range for Optimizing Systems Scheduled with Fixed Priority_. In IEEE Real-Time and Embedded Technology and Applications Symposium, pages 283–294, 2018, https://doi.org/10.1109/RTAS.2018.00036.
<br /><br />
\[3] R. Davis and A. Burns. _Improved Priority Assignment for Global Fixed Priority Pre-emptive Scheduling in Multiprocessor Real-Time Systems_. Real-Time Systems, 47(1):1–40, 2011, https://doi.org/10.1007/s11241-010-9106-5
<br /><br />
\[4] N. Guan, M. Stigge, W. Yi, and G. Yu. _New Response Time Bounds for Fixed Priority Multiprocessor Scheduling_. In 30th IEEE Real-Time Systems Symposium, pages 387–397, 2009, https://doi.org/10.1109/RTSS.2009.11.
<br /><br />
\[5] Q. Zhou, G. Li, and J. Li. _Improved Carry-in Workload Estimation for Global Multiprocessor Scheduling_. IEEE Transactions on Parallel and Distributed Systems, 28(9):2527–2538, 2017, https://doi.org/10.1109/TPDS.2017.2679195.
<br /><br />
\[6] Q. Zhou, J. Li, and G. Li. _Excluding Parallel Execution to Improve Global Fixed Priority Response Time Analysis_. ACM Transactions on Embedded Computing Systems, 20(5s):1–24, 2021, https://doi.org/10.1145/3477035
