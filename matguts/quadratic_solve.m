function [ x1,x2 ] = quadratic_solve( a,b,c )
%[ x1,x2 ] = quadratic_solve( a,b,c )
%   quick function to solve a quadratic equation given the coefficients of
%   the equation ax^2 + bx + c = 0

x1 = (-b + sqrt(b.^2 - 4*a.*c))./(2*a); 
x2 = (-b - sqrt(b.^2 - 4*a.*c))./(2*a); 

end

