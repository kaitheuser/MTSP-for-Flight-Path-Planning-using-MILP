clear all
clc

%% Flight Path Settings
num_Salesman = 5;
depot_State = 'Kansas';

%% Plot the USA map together with the visiting points
[USA_states] = USA_Map_Plot(depot_State);

%% Run the mTSP algorithm
[optimized_paths, path_seq_ij, paths_info, states_info] = mTSP(USA_states,num_Salesman, depot_State);

%% Determine each salesman route
% Initialize the route as an empty array
route = [];
% Initialize an array that stores optimized path_ij
pathOrder = zeros(length(optimized_paths), 2);
% Initialize total_path_costs
total_path_costs = zeros(1, num_Salesman);
% Get the optimized path_ij
for idx = 1:length(optimized_paths)
    pathOrder(idx, : ) = path_seq_ij(optimized_paths(idx), :);
end
% Get each salesman route
for idx = 1:num_Salesman
    % Find i with point 1
    iidx = find(pathOrder(:,1) == 1);
    fp = pathOrder(iidx(idx), 1);
    np = pathOrder(iidx(idx), 2);
    % Store salesman first path_1j
    salesman(idx).seq = [fp np optimized_paths(iidx(idx))];
    % Calculate path of the cost
    total_path_costs(1, idx) = total_path_costs(1, idx) + paths_info(salesman(idx).seq(1,3)).ij_cost(1);
    counter = 2;
    % Store the salesman path_ij until reaches the depot point
    while 1
        n_path = optimized_paths(find(pathOrder(:,1) == np));
        fp = pathOrder(find(pathOrder(:,1) == np), 1);
        np = pathOrder(find(pathOrder(:,1) == np), 2);
        salesman(idx).seq = [salesman(idx).seq; fp np n_path];
        % Calculate path of the cost
        total_path_costs(1, idx) = total_path_costs(1, idx) + paths_info(salesman(idx).seq(counter,3)).ij_cost(1);
        counter = counter + 1;
        fp = np;
        % Depot point reach at j
        if np == 1
            % Break while loop
            break
        end
    end
    % Store path coordinates
    salesman(idx).path = [paths_info(salesman(idx).seq(1,3)).i_coords(1, 1), paths_info(salesman(idx).seq(1,3)).i_coords(1, 2)];
    % Loop through all the paths
    for j = 1:size(salesman(idx).seq(:,1), 1)
        salesman(idx).path = [salesman(idx).path;
                              paths_info(salesman(idx).seq(j,3)).j_coords(1, 1), paths_info(salesman(idx).seq(j,3)).j_coords(1, 2)];
    end
end

%% Display Planned Flight Path
% Plot the USA map together with the visiting points
[USA_states] = USA_Map_Plot(depot_State);
% Color line index
lc_idx = 1;

% Plot optimized path
for idx = 1:num_Salesman
    line_colors = ['m', 'c', 'w', 'b', 'k'];
    for j = 2:size(salesman(idx).seq(:,1), 1) + 1
        plot_trajectory(salesman(idx).path(j-1, :), salesman(idx).path(j, :), line_colors(lc_idx));
    end
    lc_idx = lc_idx + 1;
    % Reset color index
    if lc_idx > 5
        lc_idx = 1;
    end
end

%% Calculate the total path costs
% U.S. UAV fuel cost based on the U.S. airline fuel cost in year 2021.
% https://www.statista.com/statistics/197689/us-airline-fuel-cost-since-2004/
fuel_price = 1.98; % USD per gallon
% https://apps.dtic.mil/sti/pdfs/ADA500380.pdf
fuel_consumption_rate = 4.2; % gallon per hour (cruise)
% UAV drone speed
% https://www.thenationalnews.com/world/2021/12/24/in-2021-we-saw-the-future-of-drone-warfare-bigger-faster-and-better-armed/#:~:text=Today%2C%20unmanned%20aircraft%20can%20fly,like%20a%20biplane%20in%20comparison.
drone_speed = 220; % km/h

%% Load the total cost of the single TSP where the starting state is California with 1 salesman.
fprintf("\nTSP Cost and Mission Completion Duration Analysis with 1 Salesman Starting at the California State\n")
fprintf("---------------------------------------------------------------------------------------------------\n\n")
total_path_len_Cali_1sm = load("mTSP_1_salesmen_California.mat");

% Calculate the mission cost in USD.
total_path_cost_Cali_1sm = total_path_len_Cali_1sm.total_path_costs(1) .* (fuel_price * fuel_consumption_rate / drone_speed);
fprintf("Total Cost: USD %3.2f \n", total_path_cost_Cali_1sm)

% Calculate the mission duration in hours.
mission_duration_Cali_1sm = total_path_len_Cali_1sm.total_path_costs(1) ./ drone_speed;
fprintf("Mission Completion Duration: %3.2f hours\n\n", mission_duration_Cali_1sm)

%% Load the total cost of the mTSP where the starting state is California with 5 salesmen.
fprintf("\nMTSP Cost and Mission Completion Duration Analysis with 5 Salesmen Starting at the California State\n")
fprintf("---------------------------------------------------------------------------------------------------\n\n")
total_path_len_Cali_5sm = load("mTSP_5_salesmen_California.mat");

% Calculate the mission cost in USD.
total_path_cost_Cali_5sm = total_path_len_Cali_5sm.total_path_costs(1,:) .* (fuel_price * fuel_consumption_rate / drone_speed);
% Calculate each salesman cost
salesman_1_cost = total_path_cost_Cali_5sm(1,1);
salesman_2_cost = total_path_cost_Cali_5sm(1,2);
salesman_3_cost = total_path_cost_Cali_5sm(1,3);
salesman_4_cost = total_path_cost_Cali_5sm(1,4);
salesman_5_cost = total_path_cost_Cali_5sm(1,5);
% Calculate total cost
total_cost_Cali_5sm = sum(total_path_cost_Cali_5sm);
fprintf("Salesman 1 Cost: USD %3.2f \n", salesman_1_cost)
fprintf("Salesman 2 Cost: USD %3.2f \n", salesman_2_cost)
fprintf("Salesman 3 Cost: USD %3.2f \n", salesman_3_cost)
fprintf("Salesman 4 Cost: USD %3.2f \n", salesman_4_cost)
fprintf("Salesman 5 Cost: USD %3.2f \n", salesman_5_cost)
fprintf("Total Cost: USD %3.2f \n\n", total_cost_Cali_5sm)

% Calculate the mission duration in hours.
mission_duration_Cali_5sm = total_path_len_Cali_5sm.total_path_costs(1,:) ./ drone_speed;
% Calculate each salesman mission duration
salesman_1_mission_duration = mission_duration_Cali_5sm(1,1);
salesman_2_mission_duration = mission_duration_Cali_5sm(1,2);
salesman_3_mission_duration = mission_duration_Cali_5sm(1,3);
salesman_4_mission_duration = mission_duration_Cali_5sm(1,4);
salesman_5_mission_duration = mission_duration_Cali_5sm(1,5);
% Select the max duration is the mission overall duration
overall_mission_duration = max(mission_duration_Cali_5sm);
fprintf("Salesman 1 Mission Completion Duration: %3.2f hours\n", salesman_1_mission_duration)
fprintf("Salesman 2 Mission Completion Duration: %3.2f hours\n", salesman_2_mission_duration)
fprintf("Salesman 3 Mission Completion Duration: %3.2f hours\n", salesman_3_mission_duration)
fprintf("Salesman 4 Mission Completion Duration: %3.2f hours\n", salesman_4_mission_duration)
fprintf("Salesman 5 Mission Completion Duration: %3.2f hours\n", salesman_5_mission_duration)
fprintf("Mission Completion Duration: %3.2f hours\n\n", overall_mission_duration)

%% Load the total cost of the single TSP where the starting state is Kansas with 1 salesman.
fprintf("\nTSP Cost and Mission Completion Duration Analysis with 1 Salesman Starting at the Kansas State\n")
fprintf("------------------------------------------------------------------------------------------------\n\n")
total_path_len_Kan_1sm = load("mTSP_1_salesmen_Kansas.mat");

% Calculate the mission cost in USD.
total_path_cost_Kan_1sm = total_path_len_Kan_1sm.total_path_costs(1) .* (fuel_price * fuel_consumption_rate / drone_speed);
fprintf("Total Cost: USD %3.2f \n", total_path_cost_Kan_1sm)

% Calculate the mission duration in hours.
mission_duration_Kan_1sm = total_path_len_Kan_1sm.total_path_costs(1) ./ drone_speed;
fprintf("Mission Completion Duration: %3.2f hours\n\n", mission_duration_Kan_1sm)

%% Load the total cost of the mTSP where the starting state is Kansas with 5 salesmen.
fprintf("\nMTSP Cost and Mission Completion Duration Analysis with 5 Salesmen Starting at the Kansas State\n")
fprintf("-----------------------------------------------------------------------------------------------\n\n")
total_path_len_Kan_5sm = load("mTSP_5_salesmen_Kansas.mat");

% Calculate the mission cost in USD.
total_path_cost_Kan_5sm = total_path_len_Kan_5sm.total_path_costs(1,:) .* (fuel_price * fuel_consumption_rate / drone_speed);
% Calculate each salesman cost
salesman_1_cost = total_path_cost_Kan_5sm(1,1);
salesman_2_cost = total_path_cost_Kan_5sm(1,2);
salesman_3_cost = total_path_cost_Kan_5sm(1,3);
salesman_4_cost = total_path_cost_Kan_5sm(1,4);
salesman_5_cost = total_path_cost_Kan_5sm(1,5);
% Calculate total cost
total_cost_Kan_5sm = sum(total_path_cost_Kan_5sm);
fprintf("Salesman 1 Cost: USD %3.2f \n", salesman_1_cost)
fprintf("Salesman 2 Cost: USD %3.2f \n", salesman_2_cost)
fprintf("Salesman 3 Cost: USD %3.2f \n", salesman_3_cost)
fprintf("Salesman 4 Cost: USD %3.2f \n", salesman_4_cost)
fprintf("Salesman 5 Cost: USD %3.2f \n", salesman_5_cost)
fprintf("Total Cost: USD %3.2f \n\n", total_cost_Kan_5sm)

% Calculate the mission duration in hours.
mission_duration_Kan_5sm = total_path_len_Kan_5sm.total_path_costs(1,:) ./ drone_speed;
% Calculate each salesman mission duration
salesman_1_mission_duration = mission_duration_Kan_5sm(1,1);
salesman_2_mission_duration = mission_duration_Kan_5sm(1,2);
salesman_3_mission_duration = mission_duration_Kan_5sm(1,3);
salesman_4_mission_duration = mission_duration_Kan_5sm(1,4);
salesman_5_mission_duration = mission_duration_Kan_5sm(1,5);
% Select the max duration is the mission overall duration
overall_mission_duration = max(mission_duration_Kan_5sm);
fprintf("Salesman 1 Mission Completion Duration: %3.2f hours\n", salesman_1_mission_duration)
fprintf("Salesman 2 Mission Completion Duration: %3.2f hours\n", salesman_2_mission_duration)
fprintf("Salesman 3 Mission Completion Duration: %3.2f hours\n", salesman_3_mission_duration)
fprintf("Salesman 4 Mission Completion Duration: %3.2f hours\n", salesman_4_mission_duration)
fprintf("Salesman 5 Mission Completion Duration: %3.2f hours\n", salesman_5_mission_duration)
fprintf("Mission Completion Duration: %3.2f hours\n\n", overall_mission_duration)
