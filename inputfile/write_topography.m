function void = write_topography(X,Y,Z,filename,my,mz)

    % A simple subroutine to write fault profiles to a file                                                   
    

    % my,mz                                        !                                    
    %file_handle                                    file handle                        
    % X                              ! the grid                          \
     
    void = 1;
    
    f = fopen([filename, '.dat'],'w');
    
    for k=1:mz
        
        for j=1:my
       
          fprintf(f,'%21.12E %21.12E %21.12E\n',X(j, k),Y(j, k),Z(j, k));
       end
    end 
    
    fclose(f);