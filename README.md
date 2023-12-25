# Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range

Code repository for Priority Assignment for Global Fixed Priority Scheduling on Multiprocessor Using Response Time Estimation Range.

## Priority Assignment Methods

## Important Note
- THe

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
- It can be moved to the **Results** directory which has sub directories with the results of other runs. 
- One of these sub directories can be duplicated, the *.txt* file in the duplicate replaced with the new one and the *.py* file can be modified to obtain the required graphs. 

## Modifying Priority Orderings (POs)

- Open ***RandomSystemsTest.cpp***.
- *TestUSweep()*:
  - The function *TestUSweep()*, specifically the part from Line 101 to Line 203, is modified to change which POs are computed.
  - The part of the function from Line 204 to the end of the file should be modified to change which POs are written into the Result.txt file for each taskset.
- *MarshallUSweepResult()*:
  - The function *MarshallUSweepResult()* in the same file is reads all the Result.txt files and compiles them into a single file. 
  - If the number of POs used has been modified in the *TestUSweep()* function, the number in the condition of the inner `for` loops should be modified.
