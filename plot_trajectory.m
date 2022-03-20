function [] = plot_trajectory(pt_i, pt_j, line_color)

% Plot a trajectory from point i (start point) to point j (destination point)
%
% Inputs : 1.) pt_i - Starting Coordinates (lat, long) 1 x 2 array
%          2.) pt_j - Destination Coordinates (lat, long) - 1 x 2 array
%          3.) line_color (string) - path color, e.g.: 'm-'
%
% Ouput  : None

% Plot Harvesine (Great Circle) Trajectory
[gclat,gclong] = track2('gc',pt_i(1), pt_i(2),...
                             pt_j(1), pt_j(2));
plotm(gclat,gclong,line_color)

end

