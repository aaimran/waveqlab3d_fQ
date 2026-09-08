function [X] = Gen_Surf_Mesh(mx,my,XL,XR,XB,XT)

   


   hq = 1/(mx-1);
   hr = 1/(my-1);

   X = zeros(my, mx);

   for j = 1:my
      for i = 1:mx
      	q = (i-1)*hq;
        r = (j-1)*hr;
        
        X(j,i) = (1-q)*XL(j)+q*XR(j)+(1-r)*XB(i)+r*XT(i)-...
               (1-q)*(1-r)*XL(1)-q*(1-r)*XR(1)-r*(1-q)*XT(1)-...
               (r*q)*XT(mx);


               
      end
      
      
   end
   
