# Multiple Traveling Salesman Problem for Flight Path Planning using Mixed-Integer Linear Programming
## Flight Path Planning Project
The mixed-integer linear programming (MILP) algorithm is implemented to solve the flight path planning problem efficiently with multiple unmanned aerial vehicles (UAVs), which can also be formulated as a multiple traveling salesman problem (mTSP).

## Flight Plan for Single TSP:
![California State As Depot](./Results/mTSP_1_Cali.png)
![Kansas State As Depot](./Results/mTSP_1_Kan.png)

## Flight Plan for mTSP:
![California State As Depot](./Results/mTSP_5_Cali.png)
![Kansas State As Depot](./Results/mTSP_5_Kan.png)


## How to Run "Visual-Inertial SLAM"?
1.) Open up the "main_Path_Planning.m" with MATLAB 2021a or newer.

2.) Run the "main_Path_Planning.m" code to perform MILP to solve mTSP or single TSP.

3.) It will output an image that has the optimized paths on the United States map.

## Parameters that Can be Changed:
1.) num_Salesman (Line 5)
- integer type
- Number of Salesmen
- Recommended input (1 - 5)

2.) depot_State (Line 6)
- string type
- The state name of the depot point
- Only input a state name that is a Contiguous US State.

## "main_Path_Planning.m" Description:
- Plot the visiting points/states and the Contiguous US States map by calling the 'USA_Map_Plot.m' function.
- Perform the MILP to solve mTSP by calling the 'mTSP.m' function.
- Determine each salesman route.
- Display the optimized flight paths by calling the 'plot_trajectory.m' and 'USA_Map_Plot.m' functions.
- Calculate the UAVs' total operating cost, and the overall mission completion duration

## "USA_Map_Plot.m" Function Description:
- Plot the Contiguous US States map.

## "mTSP.m" Function Description:
- Applied mixed-integer linear programming (MILP to solve multiple traveling salesman problem (mTSP).

## "plot_trajectory.m" Function Description:
- Plot a trajectory from point i (start point) to point j (destination point)

## "haversine_distance.m" Function Description:
- Calculate the Haversine distance between 2 coordinates in km.
