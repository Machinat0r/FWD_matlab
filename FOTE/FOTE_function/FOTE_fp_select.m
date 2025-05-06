function [h,focus_points] = FOTE_fp_select(X_coords, Y_coords, Z_coords, coord_type)

if nargin == 3
    if length(X_coords) == length(Y_coords) == length(Y_coords)
       focus_points = zeros(length(X_coords),3);
       focus_points(:,1) = X_coords;
       focus_points(:,2) = Y_coords;
       focus_points(:,3) = Z_coords;
    else
        coord_type = 'Cartesian';
    end
end
m = 1;
if nargin == 4
    focus_points = zeros(length(X_coords)*length(Y_coords)*length(Z_coords),3);
    
    switch coord_type
        case 'Cartesian'
            for i = 1:length(X_coords)
                for j = 1:length(Y_coords)
                    for k = 1:length(Z_coords)
                        focus_points(m,1) = X_coords(i);
                        focus_points(m,2) = Y_coords(j);
                        focus_points(m,3) = Z_coords(k);
                        m = m+1;
                    end
                end
            end
            
        case 'Cylindrical'
            for i = 1:length(X_coords)
                for j = 1:length(Y_coords)
                    for k = 1:length(Z_coords)
                        focus_points(m,1) = X_coords(i)*cosd(Y_coords(j));
                        focus_points(m,2) = X_coords(i)*sind(Y_coords(j));
                        focus_points(m,3) = Z_coords(k);
                        m = m+1;
                    end
                end
            end
        case 'Spherical'
            for i = 1:length(X_coords)
                for j = 1:length(Y_coords)
                    for k = 1:length(Z_coords)
                        focus_points(m,1) = X_coords(i)*sind(Y_coords(j))*cosd(Z_coords(k));
                        focus_points(m,2) = X_coords(i)*sind(Y_coords(j))*sind(Z_coords(k));
                        focus_points(m,3) = X_coords(i)*cosd(Y_coords(j));
                        m = m+1;
                    end
                end
            end
    end
end

if nargin == 2
    d1 = size(X_coords);
    d2 = size(Y_coords);
    if d1(2) == d2(2) == 3
        focus_points = [X_coords; Y_coords];
    else
        disp('Please check the formats of arguments!')
        edit FOTE_fp_select
        return
    end   
end

h = [];
% h(1) = axes('position',[0.1 0.1 0.8 0.8]);

% h = plot3(gca, focus_points(:,1),focus_points(:,2),focus_points(:,3), 'bo', 'Linewidth',1, ...
% 'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);


end