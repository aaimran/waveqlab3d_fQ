function [gg] = interpolate_surface(X1, Y1, Z1, y, z, my, mz)

    

    
    
    stenc = 1;
    ndp = 3;
   
                
    for kk = 1:mz
       for jj = 1:my
           
          dy = Y1(2, 1) - Y1(1,1);
          dz = Z1(1, 2) - Z1(1,1);
          
          if (jj > 1)
             dy = Y1(jj, kk) - Y1(jj-1,kk); 
          end
          
          if (kk > 1)
             dz = Z1(jj, kk) - Z1(jj,kk-1); 
          end
              

          r0 = sqrt(dy^2 + dz^2);

          r = sqrt((y-Y1(jj,kk))^2 + (z-Z1(jj,kk))^2);

          if (r <= 2.0*r0)
             
             if (abs(y-Y1(jj,kk)) <= dy)
                
                nmy = jj-stenc;
                npy = jj+stenc;
                
                
                if (jj == 1) 
                   
                   nmy = jj;
                   npy = jj+ndp-1;
                   
                end

                if (jj == 2)
                   
                   nmy = jj-1;
                   npy = jj+ndp-2;
                   
                end 

                 if (jj == 3)
                   
                   nmy = jj-2;
                   npy = jj+ndp-3;
                   
                end 

                if (jj == my-2)
                   
                   nmy = jj-(ndp-3);
                   npy = jj+2;
                   
                end
                
                if (jj == my-1)
                   
                   nmy = jj-(ndp-2);
                   npy = jj+1;
                   
                end

                if (jj == my);
                   
                   nmy = jj-(ndp-1);
                   npy = jj;
                   
                end;
                
                
                if  (abs(z-Z1(jj,kk)) <= dz)
                   
                   nmz = kk-stenc;
                   npz = kk+stenc;
                   
                   if (kk == 1)
                      
                      nmz = kk;
                      npz = kk+ndp-1;
                      
                   end

                   if (kk == 2)
                      
                      nmz = kk-1;
                      npz = kk+ndp-2;
                      
                   end

                   if (kk == 3)
                      
                      nmz = kk-2;
                      npz = kk+ndp-3;
                      
                   end
                   
                   if (kk == mz)
                      
                      nmz = kk-(ndp-1);
                      npz = kk;
                      
                   end

                   if (kk == mz-1)
                      
                      nmz = kk-(ndp-2);
                      npz = kk+1;
                      
                   end

                   if (kk == mz-2)
                      
                      nmz = kk-(ndp-3);
                      npz = kk+2;
                      
                   end
                     
                   
                   gg = interpol2d_dG(nmy, npy, nmz, npz, y, z, Y1(:,kk), Z1(jj, :), X1(:, :));
                                    
               
                end
                %==============================================
             end
          end

           %print *,  jj, kk, Y1(jj,kk), Z1(jj, kk)
       end
    end
    
end
