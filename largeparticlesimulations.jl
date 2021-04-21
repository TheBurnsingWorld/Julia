
using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Emergent1_6_Angle_45.0_NP_1000_NV_4_Rho_0.100_Config.dat")
#h=readdlm("Emergent1_6_Angle_45.0_NP_1000_NV_4_Rho_0.100_LargestCluster.dat")
f=f[78:end,1:9]
#timescale=(f[end,end]/4)+1
#timescale=trunc(Int64,timescale)
Two=Array{Int32}(undef, 1592)
Three=Array{Int32}(undef, 1592)
n=1003
for i in 1:1592
    Two[i]=(n-1)+((i-1)*(n))
    Three[i]=Two[i]+1
end
#find the positions of the large particles
#get the values from column 1,2
TwoPos=Array{Float32}(undef,1592,2)
ThreePos=Array{Float32}(undef,1592,2)
Distance=Array{Float32}(undef,1592)
Difference=Array{Float32}(undef,1592,2)
Length=200
MinDistance=7.07
for i in 1:1592
    TwoPos[i,1]=f[Two[i],1]
    ThreePos[i,1]=f[Three[i],1]
    TwoPos[i,2]=f[Two[i],2]
    ThreePos[i,2]=f[Three[i],2]


rminmat=zeros(3,3)
for p in [-1,0,1]
    for q in [-1,0,1]
        function g(p,q)

            (sqrt((((ThreePos[i,1]-TwoPos[i,1]+(p*Length))^(2))+(((ThreePos[i,2]-TwoPos[i,2]+(q*Length))^(2)))))) #CHANGE THIS LINE SO THAT THE DISTANCE IS NOT LESS THAN IT SHOULD BE FOR THE LARGE PARTICLES
        end
        if g(p,q)>MinDistance
            rminmat[(p+2),(q+2)]=g(p,q) # we want to find the smallest value that is larger than the smallest acceptable distance
        else
            rminmat[(p+2),(q+2)]=g(p,q)+100000
        end
        # define the global smallest acceptable distance.
        # find the
        #if g(p,q)<
    end
end
rmin=minimum(rminmat)
vals=argmin(rminmat)
p1=vals[1]
q1=vals[2]
Difference[i,:]=[(ThreePos[i,1]-TwoPos[i,1]+((p1-2)*Length)),(ThreePos[i,2]-TwoPos[i,2]+((q1-2)*Length))]

end

#Difference=TwoPos-ThreePos
for i in 1:1592
    Distance[i]=norm(Difference[i,:])
end

#finding the histogram
N=zeros(1592)
Distance1=Array{Int32}(undef,1592)
for i in 1:1592
    Distance1[i]=Int32(ceil(Distance[i]))
end
for i in 1:1592
    N[Distance1[i]]=N[Distance1[i]]+1
end
N1=zeros(1592)
for i in 1:1592
    N1[i]=N[i]/i
end
