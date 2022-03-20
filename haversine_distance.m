function [distance] = haversine_distance(Lat_1, Lon_1, Lat_2, Lon_2)
% Calculate the distance between 2 coordinates in km.
%
% Inputs : 1.) Lat_1 - point i lattitude
%          2.) Lon_1 - point i longitude
%          3.) Lat_2 - point j lattitude
%          4.) Lon_2 - point j longitude
%
% Ouputs : 1.) distance - harvesine distance

% Lattitude and Longitude difference
d_Lat = (Lat_2 - Lat_1) * pi / 180;
d_Lon = (Lon_2 - Lon_1) * pi / 180;

% Convert to radian
Lat_1 = Lat_1 * pi / 180;
Lat_2 = Lat_2 * pi / 180;

% Calculate the coefficient a
a = (sin(d_Lat/2))^2 + cos(Lat_1) * cos(Lat_2) * (sin(d_Lon/2))^2;
% Calculate the coefficient c
c = 2 * asin(sqrt(a));
% Radius of the earth [km]
r = 6371;

% Calculate the distance
distance = c * r;

end

