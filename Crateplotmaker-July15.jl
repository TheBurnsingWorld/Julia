using Plots
using DelimitedFiles
cd("E:/Research Data/Computation/Active/Calculate_Rate_Constants April 2019")
f=readdlm("Test_May_15_TC1000_0_matrices.dat")
g=f[11:end,1:3]
a=zeros(1000,1000)
r=g[:,1]
c=g[:,2]
v=g[:,3]

function h(x,y)

              (1000)*(y-1)+x

              end

for i in 1:1000

                for j in 1:1000

                    a[i,j]=g[ h(r[i],r[j]),3]

                    end

                    end


b=zeros(1000,1000)

for i in 1:1000
       for j in 1:1000
       if a[i,j]>0
       b[i,j]=log(a[i,j])
       end
       end
       end      
