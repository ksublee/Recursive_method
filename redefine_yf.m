function [y, dy, f] = redefine_yf(x, y, f, Num_y, method)
    dy = y(2) - y(1);
    if length(y) > Num_y+ 1
        dy = (max(y)-min(y))/Num_y;
        new_y = min(y):dy:max(y);
        [mesh_y, mesh_x] = meshgrid(new_y, x);
        f = interp2(y, x, f, mesh_y, mesh_x, method);
        y = new_y;
    end
end