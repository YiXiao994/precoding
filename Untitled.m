       for m = 1:settings.num_of_Antenna
          for n = 1:settings.num_of_Antenna
              if m == n
                  f_X(m,n) = 2 * X(m,n) - sum(X(m,:));
              else
                  f_X(m,n) = X(m,n);
              end
              
          end           
       end