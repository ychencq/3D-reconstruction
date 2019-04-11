clear
I = imread('test2.png');
 I = double(rgb2gray(I))/255;
I = double(I)/255;

[numRow,numCloumn] = size(I);


% for we need to know the 2 deflection angles in an actual problem
% here we set £ps(phiS) as -pi/2 (-90) and £cs(thetaS) as pi/2 (90) 

phiS = -(pi/2);
thetaS = (pi/2);

% and conversion matrix   :       R [cos£ps*cos£cs cos£ps*sin£cs -sin£ps]
%                                   [-sin£cs         cos£cs      0   ]
%                                   [cos£ps*cos£cs sin£ps*sin£cs  cos£ps]

% hecence we know that the R here should be:
%R = [0 0 1; -1 0 0; 0 -1 0];

%if we product R with the greyscale image we can get a true grayscale
% image E      (E = [x;y;z]*R)


% and the function of E is to elimate the effect of deflection angles of
% the light, let we consider the situation that we are looking the photo
% orthophorially

% now we want to get the gradient of E¡A and it is exactly:

E3d = zeros(numRow,numCloumn,1000);
% notice that for the  y-axis & z-axis the poistive is actually negative number!!!!!
for i = 1:numRow
    for j = 1:numCloumn
    E3d(1,i,j) = I(i,j);
    end
end



Izx = zeros(numRow,numCloumn);
Izy = zeros(numRow,numCloumn);

for i = 1:numRow
    for j = 2:numCloumn
    Izx(i,j) = I(i,j)-I(i,j-1);
    end
end

for i = 1:numRow
    Izx(i,1) = 0;
end

for n = 1:numCloumn
    for m = 2:numRow
        Izy(m,n) = I(m,n)-I(m-1,n);
    end
end

for n = 1:numCloumn
    Izy(1,n) = 0;
end

% and we can transfer Izx Izy to Ezx Ezy
 Ezx = Izx*cos(phiS)*cos(thetaS) + Izy*cos(phiS)*sin(thetaS);
 Ezy = -Izx*sin(thetaS) + Izy*cos(thetaS);
 Emax = max(E3d(:));
% for we got the gradient of the E¡A now we try to get the dip angle £p (Phi)
% and  drift angle £c ¡]Theta¡^

% the calculation is based on Lambert Cosin Law:
%           Ei = Emax*cos£pi => £pi=arccos(Ei/Emax)   and       £ci = arctan(Ezy/Ezx);
%still notice the actual positions these values should be 
phi =   acos(E3d/Emax);
theta = atan(Ezy./Ezx);
% here notice that thetai here is 2 2Dmatrix for xy can determine it.
theta3D = zeros(numRow,numCloumn,1000);
for i = 1:numRow
    for j = 1:numCloumn
        for k = 1:1000
            theta3D(i,j,k) = theta(i,j);
        end
    end
end

%  to get the normal vector (nx,ny,nz) and the height
reflection = 8.35;
%  nx = sin(phi).*cos(theta3D);
%  ny = sin(phi).*sin(theta3D);
 nz = cos(phi);
 height = nz*reflection;
 
 %now we transfer back to the orginial coordinate system using the normal vector
 %with R' = [-cos£ps*cos£cs    sin£ps     sin£ps*cos£cs]
 %          [ cos£ps*sin£cs    cos£cs    -sin£ps*sin£cs]
 %          [    sin£ps        0           cos£ps   ]
 
 % hence R' = [0 -1 0; 0 0 1; -1 0 0]
 % still remeber the index problem before!!!
 
 Inew = zeros(numRow,numCloumn,1000);
 
 for i = 1:numRow
     for j = 1:numCloumn
         for k = 1:1000
         Inew(i,j,k) =  height(j,k,i);
         end
     end
 end
 
 Inew = Inew(:,:,1);
 
 
 
[x,y] = size(Inew);
[X,Y] = meshgrid(1:y,1:x);

s =  surf(X,Y,Inew,'FaceAlpha',0.5)
s.EdgeColor = 'none';
colorbar
 
 
 



