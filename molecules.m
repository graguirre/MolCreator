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
	bond = 1.4; % C-C
	a = bond/(2*tan( pi/p ));
	r = bond/(2*sin( pi/p ));
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
function print(v)
	for i=1:length(v)
		printf("C %.6f %.6f 0.0\n ", v(1,i), v(2,i))
	end
end
% create pentagon
[a5,r5]=apotheme(5);
x5=polygon(5,r5);
%create heptagon
[a7,r7]=apotheme(7);
x7=polygon(7,r7);
% pentagon translation matrix
x5 = translate(x5, a5+a7, 0);
x5 = rotate(x5, pi/7);
% attach azulene molecule
C10H8 = [x7 x5]; % there are 2 repeated atoms
C10H8(:,[1 7]) = []; % erase repeated atoms

%hold on; % uncomment if debugging
C10H8 = translate(C10H8, -C10H8(1,3), -C10H8(2,3));
%plot(C10H8(1,:),C10H8(2,:),'r.'); % uncomment if debuigging

print(C10H8);
% get rotation
v = [C10H8(1,10)-C10H8(1,9) C10H8(2,10)-C10H8(2,9)];
alpha = atan( v(1)/v(2) );

t=[C10H8(1,9) C10H8(2,9)];

C10H8(:,[2 3]) = []; % erase repeated atoms
for i=1:32
	C10H8 = rotate(C10H8,-alpha);
	C10H8 = translate(C10H8, t(1), t(2));
	%plot(C10H8(1,:),C10H8(2,:),'b.'); % uncomment if debugging
	print(C10H8);
end

%input("Press enter to continue...");
