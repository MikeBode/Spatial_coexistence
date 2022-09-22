% Mathematical Modelling in Ecology Assessment 2 Code
% By Dylan Maher

function Shape_Index = Code_Dylan_Maher(IslandOutline)


%%% Calculate the perimeter of the island in kilometres. Note that:
% - Longitude is given by IslandOutline(:,1)
% - Latitude is given by IslandOutline(:,2)
% - The following method will be used to compute the distance between two points on the surface of the earth:
%   dlon = longitude2 - longitude1
%   dlat = latitude2 - latitude1
%   b = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
%   c = 2*atan(sqrt(b), sqrt(1-b))
%   distance_between_two_points = R*c

R = 6371; % The radius of the earth in kilometres
points = size(IslandOutline, 1); % The number of points on the edge of the island that are given

% Find the distance between the first and last points to initialise the perimeter 
dlon = IslandOutline(1,1) - IslandOutline(points,1);
dlat = IslandOutline(1,2) - IslandOutline(points,2);
b = (sind(dlat/2)^2) + (cosd(IslandOutline(points,2))*cosd(IslandOutline(1,2))*(sind(dlon/2)^2));
c = 2*atan2d(sqrt(b),sqrt(1-b));
Perimeter = R*c;

% Calculate the rest of the perimeter (add the distance between the first and second points, add the distance between the second and third
% points,..., add the distance between the second last and last points)
for i=1:(points-1)
    dlon = IslandOutline(i+1,1) - IslandOutline(i,1);
    dlat = IslandOutline(i+1,2) - IslandOutline(i,2);
    b = (sind(dlat/2)^2) + (cosd(IslandOutline(i,2))*cosd(IslandOutline(i+1,2))*(sind(dlon/2)^2));
    c = 2*atan2d(sqrt(b),sqrt(1-b));
    Perimeter = Perimeter + (R*c);
end


%%% Calculate the area of the shape of the island in kilometres
Area = areaint(IslandOutline(:,2), IslandOutline(:,1), earthRadius('km'));


%%% Calculate the shape index for the island
Shape_Index = Perimeter/(2*sqrt(pi*Area));