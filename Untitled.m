clear all;
tic
t1 = cputime
i = 0
for a = 1:100000
   i = i + 1; 
end
t2 = cputime - t1
toc
