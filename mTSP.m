function [optimized_paths, path_seq_ij, paths_info, states_info] = mTSP(USA_states, num_Salesman, depot_State)
% Applied mixed-integer linear programming (MILP to solve multiple
% traveling salesman problem (mTSP).
%
% Inputs : 1.) USA_states - an object that contains states information.
%          2.) num_Salesman - number of salesmen to complete the mission.
%          3.) depot_State - initial / last end point of the mission.
%
% Ouputs : 1.) optimized_paths - optimized flight paths
%          2.) path_seq_ij - all possible path_ij in an array form
%          3.) paths_info - all possible paths information 
%          4.) states_info - visiting states information.

% Get state names, lat, and lon
state_Names = {USA_states.Name};
state_Lat = {USA_states.LabelLat};
state_Lon = {USA_states.LabelLon};

% Get number of visiting points
num_states = length(state_Names);

%% Initialize an object that stores stateID, state name, and coordinates called 'states_info'
% Identify the depot state and store it as the first visiting point
for state = 1 : num_states
    % If depot state is found
    if strlength(state_Names{state}) == strlength(depot_State)
        if state_Names{state} == depot_State
            % Store information
            states_info(1).pt_id = 1;
            states_info(1).name = state_Names{state};
            states_info(1).lat = state_Lat{state};
            states_info(1).lon = state_Lon{state};
            % Break for loop
            break
        end 
    end
end

% Starts the state ID with 2 since 1 is taken by the depot state
ID_count = 2;

% Initialize other visiting states
for state = 1 : num_states
    
    % If depot state is found
    if strlength(state_Names{state}) == strlength(depot_State)
        if state_Names{state} == depot_State
            % Do nothing.
            continue  
        end
    end
    
    % Store information
    states_info(ID_count).pt_id = ID_count;
    states_info(ID_count).name = state_Names{state};
    states_info(ID_count).lat = state_Lat{state};
    states_info(ID_count).lon = state_Lon{state};
    
    % Increment the count
    ID_count = ID_count + 1;
    
end

%% Store path information including path_ij, path_cost, i_coordinates, j_coordinates, i_state, j_state
% Initialize field number 
field_num = 0;

% Intialize path_seq_ij
path_seq_ij = zeros(num_states * (num_states - 1), 2);

% An array that determines and stores the minimum cost for each pair of checkpoints with its ID.
minCostIDnPath = zeros(num_states * (num_states - 1), 2);

% Starts store path information
for i = 1:num_states
    for j = 1:num_states
        % If same states
        if i == j
            % Skip this j for loop
            continue
        % Else, store the path information
        else
            % Increment the field number
            field_num = field_num + 1;
            % Initialize path_ij order
            path_order(field_num, 1) = i;
            path_order(field_num, 2) = j;
            % Intialize paths_info
            % Path ij
            paths_info(field_num).ij_seq = path_order(field_num, :);
            path_seq_ij(field_num, :) = path_order(field_num, :);
            % Path cost [km] 
            path_cost = haversine_distance(states_info(i).lat, states_info(i).lon,...
                                           states_info(j).lat, states_info(j).lon);
            paths_info(field_num).ij_cost = path_cost;
            minCostIDnPath(field_num, 1) = field_num; % Store path ID
            minCostIDnPath(field_num, 2) = path_cost; % Store path cost
            % i_coordinates
            paths_info(field_num).i_coords = [states_info(i).lat, states_info(i).lon];
            % j_coordinates
            paths_info(field_num).j_coords = [states_info(j).lat, states_info(j).lon];
            % i_state
            paths_info(field_num).i_state = states_info(i).name;
            % j_state
            paths_info(field_num).j_state = states_info(j).name;
        end
    end
end

%% Extract path costs
path_costs = cell2mat({paths_info.ij_cost})';
% Get number of possible path_ij
num_paths_ij = length(path_costs);
% path_costs size should be ((num_ij + num_states - 1) x 1)
path_costs = [path_costs; zeros(num_states - 1, 1)];

%% MATLAB TSP: Solver-Based Methods
%----------------------------------
% Equality Constraints
%*********************
% A_eq * x = b_eq 
% Initialize A matrix for equality constraints
A_eq = zeros(2 * num_states, length(num_paths_ij) + num_states - 1);
% Initialize b vector for equality constraints
b_eq = zeros(2 * num_states, 1);
    
% The first constraint enforces that all visiting points/states must be visited once.
% For departure
for stateID = 1: num_states
    
    % Find the paths that include a specific visiting point.
    % "whichPath" is a logical array ( 1 or 0 only).
    whichPath = (path_seq_ij(:,1) == stateID);

    % Include in the constraint matrix.
    A_eq(stateID, 1:num_paths_ij) = whichPath';
    A_eq(stateID, num_paths_ij + 1 : num_paths_ij + num_states - 1) = zeros(1, num_states - 1);
    
end

% "b_eq" is a (2 x number of states) x 1 matrix because there are 2 equality constraints.
% b_eq vector is the number of paths.
for stateID = 1 : num_states
    
    if stateID == 1
        
        % Equality constraint at depot is the number of salesman
        b_eq(stateID) = num_Salesman;
        
    else

        % Equality constraint at other states is 1
        b_eq(stateID) = 1;
        
    end
end

% For arrival
for stateID = 1: num_states
    
    % Find the paths/trees that include a specific checkpoint.
    % "whichPath" is a logical array ( 1 or 0 only).
    whichPath = (path_seq_ij(:,2) == stateID);
    A_eq(num_states + stateID, 1 : num_paths_ij) = whichPath';
    A_eq(num_states + stateID, num_paths_ij + 1 : num_paths_ij + num_states - 1) = zeros(1,num_states - 1); % For u_i
    
end

% "b_eq" is a (2 x number of states) x 1 matrix because there are 2 equality constraints.
% b_eq vector is the number of paths.
for stateID = 1 : num_states
    
    if stateID == 1
        
        % Equality constraint at depot is the number of salesman
        b_eq(num_states + stateID) = num_Salesman;
        
    else
        
        % Equality constraint at other states is 1
        b_eq(num_states + stateID) = 1;
        
    end
end

% Inequality Constraints
%************************
% A * x = b
% Initialize A matrix for inequality constraints
A = zeros((num_states - 1) * (num_states - 2), num_paths_ij + num_states - 1);
% Initialize b vector for inequality constraints
b = zeros((num_states - 1) * (num_states - 2), 1);

% Inequality counter
Ineqnum = 0;
% Sub-tour constraints
for i = 2 : num_states
    for j = 2 : num_states
       % Only assign when i is not equal to j
       if i ~= j
          Ineqnum =  Ineqnum + 1;
          A(Ineqnum, 1 : num_paths_ij + num_states - 1) = zeros(1, num_paths_ij + num_states - 1);
          A(Ineqnum, num_paths_ij + i - 1) = 1; % u_i
          A(Ineqnum, num_paths_ij + j - 1) = -1; % u_j
          % If i is bigger than j
          if i<j
              A(Ineqnum, (num_states - 1) * (i - 1) + j - 1) = num_states - num_Salesman + 1; % xij
          else
              A(Ineqnum,(num_states - 1) * (i - 1) + j) = num_states - num_Salesman + 1; % xij
          end
          b(Ineqnum) = num_states - num_Salesman;
       end
    end    
end

% Binary Bounds for the x_mtsp
%*****************************
% Number of decision variables (x_ij) or index of x_mtsp.
x_mtsp_index = 1 : (num_paths_ij + num_states - 1);
% Lower bound for the x_ij is zero.
lower_Bound = zeros(num_paths_ij + num_states - 1, 1);
% Lower bound for the u_i and u_j is 1.
lower_Bound(num_paths_ij + 1 : end) = ones(1, num_states - 1);
% Upper bound for the x_ij is one.
upper_Bound = ones(num_paths_ij + num_states - 1, 1);
% Upper bound for the u_i and u_j is the maximum number of destination
% points
upper_Bound(num_paths_ij + 1 : end) = num_states * ones(1, num_states - 1);

% Optimizing the solution of x_mtsp with MILP using "intlinprog"
%***************************************************************
MILP_Settings = ...
    optimoptions('intlinprog','Display','iter','Heuristics','advanced');

[x_mtsp] = ...
    intlinprog(path_costs, x_mtsp_index, A, b, A_eq, b_eq,...
               lower_Bound, upper_Bound, MILP_Settings);
           
% Round up to nearest integer        
x_mtsp = round(x_mtsp);

% Find out which x_mtsp indices are selected by the MILP
minCostIDnPathIndex = find(x_mtsp == 1);
for idx = 1:length(minCostIDnPathIndex)
    % If number of possible paths exceeded
    if minCostIDnPathIndex(idx) > num_paths_ij
        % Break the for loop
        break
    end
    % Get the index with minimum path cost
    realMinCostIDnPathIndex(idx) = minCostIDnPathIndex(idx);
end

optimized_paths = minCostIDnPath(realMinCostIDnPathIndex', 1);


end

