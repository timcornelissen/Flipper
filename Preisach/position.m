function  Y = position(X,Grids,In)
% POSITION gives the position and unit dipole vector of a specific dipole in Cartesian coordinates
% Y = POSITION(X,Grids,In) 

%% Create Y and set base x,y,z position
Y = zeros(3,2);
Y(:,1) = [X(1)*In.a+X(2)*In.a/2, X(2)*sqrt(3)/2*In.a, X(3)*In.c];

%% Set alpha for the 3 different dipoles
switch X(4)
    case 1
        alpha = Grids.Alpha(X(1),X(2),X(3));
    case 2
        alpha = Grids.Alpha(X(1),X(2),X(3))+2*pi/3;
    case 3
        alpha = Grids.Alpha(X(1),X(2),X(3))-2*pi/3;
end

%% Set dipole vector
if Grids.chirality(X(1),X(2),X(3)) %lefthanded
    Y(:,2) = [cos(alpha)*cos(In.beta); -sin(alpha)*cos(In.beta); sin(In.beta)];
else %righthanded
    Y(:,2) = [-cos(alpha)*cos(In.beta); sin(alpha)*cos(In.beta); sin(In.beta)]; 
end
    
%% Set correct x,y position
Y(1,1) = Y(1,1)+In.L*sin(alpha);%x position
Y(2,1) = Y(2,1)+In.L*cos(alpha);%y position

end
