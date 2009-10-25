C      This function just returns a double.  
C      This function serves to test if the return value
C      can be correctly captured in C.
C
       double precision function dfunction()
       dfunction = 9.87 
       return 
       end
       
       double complex function dzfunction(d,i)
       double complex temp
       double precision t(2), d, i
       equivalence (temp,t(1))
       t(1) = d
       t(2) = i
       dzfunction = temp
       return
       end
