function [USA_states] = USA_Map_Plot(start_State)

% Function that plot USA map
% Input : start_State
% Input Type : string
% Output: USA_state
% Output Type: Object

% Load the USA states shape file to obtain 
% [Field, Geometry, Bounding Box, Lon, Lat, Name, LabelLat, LabelLong]
USA_states = shaperead('usastatehi.shp','UseGeoCoords',true);

% Removes the non-Contiguous US States which are Alaska and Hawaii
USA_states(2) = []; % Alaska
USA_states(10) = []; % Hawaii

% Initialize figure 
figure('Name','Contiguous United States Map','NumberTitle','off');
ax = usamap('conus');

% Ocean Color in normalized RGB format
% ocean_Blue = [0.3010 0.7450 0.9330];
ocean_Blue = [0 0.4470 0.7410];

% Land Color in normalize RGB format
%land_Orange = [0.9290 0.6940 0.1250];
land_Green = [0.4660 0.6740 0.1880];

% Plot the Contiguous United States Map
setm(ax,'FFaceColor',ocean_Blue)
geoshow(USA_states,'FaceColor',land_Green)
title({'Contiguous United States Map', ...
    '(Multiple Traveling Salesman Problem)'})

% Get Lat and Lon of all cities
state_Lats = [USA_states.LabelLat];
state_Lons = [USA_states.LabelLon];

% Get state names
state_Names = {USA_states.Name};

% Plot Capital Cities of Each States
for city = 1 : length(state_Names)
    
    % If it is the starting state, plot different point style to differentiate with
    % other visiting states
    if strlength(state_Names{city}) == strlength(start_State)
        if state_Names{city} == start_State
            geoshow(state_Lats(city), state_Lons(city), 'DisplayType', 'Point', 'Marker', '*', 'Color', 'red');
            continue
        end
    end
    
    % Plot standard point style for other states
    geoshow(state_Lats(city), state_Lons(city), 'DisplayType', 'Point', 'Marker', 'o', 'Color', 'red');
    
end

end

