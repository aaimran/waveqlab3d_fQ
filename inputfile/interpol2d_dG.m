function [g] =  interpol2d_dG(nmx, npx, nmy, npy, x, y, xi, yj, f)
    
    %contruct the langrange interpolation g(x,y) of  f(x,y)
    
    % g(x,y): lagrange interpolant
    % nmx, nmy: lower bound
    % npx, npy: upper bound
    % xi,yj: data points
    % x,y : interpolation point
       
    g = 0.0;
    
    for j = nmy:npy
       
       for i = nmx:npx
          
          a_i = lagrange(nmx, npx, i, x, xi);
          a_j = lagrange(nmy, npy, j, y, yj);
          
          g = g + a_i*a_j*f(i,j);
          
       end
    end
    
end


  function [h] = lagrange(m, p, i, x, xi) 
    
    %Function to calculate  Lagrange polynomial for order N and polynomial
    %[nm, np] at location x.
    
    %nm: lower bound
    %np: upper bound
    %xi: data points
    %x : interpolation point
   
    
    h = 1.0;
    for j = m:p
       if (j ~= i)
          
          num = x - xi(j);
          den = xi(i) - xi(j);
          h = h*num/den;
          
       end
    end
    
  end 
  

