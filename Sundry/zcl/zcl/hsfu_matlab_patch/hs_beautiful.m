function J = hs_beautiful(m)


if nargin < 1
   m=51;
end
m=51;

%forth turn
Ncolor=m;
step=255/Ncolor;

R=[zeros(1,Ncolor)+255, 255-([1:Ncolor].*step), zeros(1,Ncolor),     zeros(1,Ncolor),        ([1:Ncolor].*step)  ];
G=[([1:Ncolor].*step),  zeros(1,Ncolor)+255,    zeros(1,Ncolor)+255, 255-([1:Ncolor].*step), zeros(1,Ncolor)     ];
B=[zeros(1,Ncolor),     zeros(1,Ncolor),        ([1:Ncolor].*step),  zeros(1,Ncolor)+255,    zeros(1,Ncolor)+255 ];

R=wrev(R);
G=wrev(G);
B=wrev(B);


RR=R./255;
GG=G./255;
BB=B./255;


J = zeros(5*m,3);
J(:,1) = RR;
J(:,2) = GG;
J(:,3) = BB;
