%
% Azelene molecule creator. (beta)
% versoin: 0.1
%
% Gonzalo Aguirre <graguirre@gmail.com>
%

% bond distance (Angstroms)
%
% C-C 1.20 -- 1.54
% C=C 1.33 
% C-H 1.06 -- 1.12
% 

% return set of coordinates centered at origin
function x=polygon(p,r)
	for i=[1:p]
		x(1,i) = r * cos( 2*pi/p*i );
		x(2,i) = r * sin( 2*pi/p*i );
	end
end
% return apotheme of a 1 size side polygon (P=#-of-sized)
function [a,r]=apotheme(p)
	bondCC = 1.35; % C-C
	a = bondCC/(2*tan( pi/p ));
	r = bondCC/(2*sin( pi/p ));
end
function x=translate(X,tx,ty)
	% create translation matrix
	T = eye(3);
	T(1,3) = tx;
	T(2,3) = ty;
	X(3,:) = ones; % append to make homnogenous transform
	x = T * X;
	x = x(1:2,:);
end
function x=rotate(X,theta)
	R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
	x = R * X;
end

%
% print atoms coordinate 
% parameters: v vector, r repeated molecule (0 first, 1 other case)
%
function print(v,r) 
	if (r == 0)
		r = 11;
	else
		r = 9;
	end
	for i=1:length(v)
		if (i < r)
			printf(" C\t%.1f\t%.6f\t%.6f\t%.6f\n", 6, v(1,i), v(2,i), 0)
		else
			printf(" H\t%.1f\t%.6f\t%.6f\t%.6f\n", 1, v(1,i), v(2,i), 0)
		end
	end
end

%
%%% Main program 
%

bondCH = 1.1;
hydro = 0;	% hydrogenate
plt = 1;	% plot
mol = 2;	% number of molecules

% create regular pentagon
[a5,r5] = apotheme(5);
C5 = polygon(5, r5);
H5 = [];	
if ( hydro == 1 )
	H5 = polygon(5, r5 + bondCH);  % hydrogenate
	H5(:,[2:5]) = [];
end

%create regular heptagon
[a7,r7] = apotheme(7);
C7 = polygon(7, r7);
H7 = [];
if ( hydro == 1)
	H7 = polygon(7, r7 + bondCH); % hydrogenate
	H7(:,[1 3 4 7]) = [];
end

% pentagon translation matrix
C5 = translate([C5 H5], a5+a7, 0);
C5 = rotate(C5, pi/7);
% attach azulene molecule
C10H8 = [C7 C5 H7]; % there are 2 repeated atoms
C10H8(:,[1 7]) = []; % erase repeated atoms

C10H8 = translate(C10H8, -C10H8(1,3), -C10H8(2,3));
if ( plt == 1 )
	hold on; % uncomment if debugging
 	plot(C10H8(1,:),C10H8(2,:),'r.'); % uncomment if debuigging
end

print(C10H8, 0);
% get rotation
v = [C10H8(1,10)-C10H8(1,9) C10H8(2,10)-C10H8(2,9)];
alpha = atan( v(1)/v(2) );

t=[C10H8(1,9) C10H8(2,9)];

C10H8(:,[2 3]) = []; % erase repeated atoms
for i=1:mol-1
	C10H8 = rotate(C10H8,-alpha);
	C10H8 = translate(C10H8, t(1), t(2));
	if ( plt == 1 )
		plot(C10H8(1,:),C10H8(2,:),'b.'); % uncomment if debugging
	end
	print(C10H8,1);
end

if ( plt == 1 )
	input("Press enter to continue...");
end
