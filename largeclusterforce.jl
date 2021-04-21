using DelimitedFiles
using Plots
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_Config.dat")
h=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_LargestCluster.dat")
#goals:
#define a set (or more than one set) of large clusters by:
#go by time steps, if largest cluster is as big as predefined limit to look at large clusters
#select that time step
# take two largest clusters, and find center of mass
#CM=(1/N)<sum x,sum y>
g=f[38:end,1:9]
limit=100
rho=f[17,4] #this should remain constant, check if not.
#define the timesteps with large clusters, or arrays from
lc=[]
#define number of timesteps,ts

#define number of particles
n=2000
j=h[37:end,1:8] # cutoffs for largest cluster data file
j1=j[:,1:2]
s=size(j1)

for i in 1:s[1]
    if j1[i,2]>=limit
        push!(lc,j1[i,1])
    end
end
#remove timestamps from that if there are less than 10 continuous datapoints
# if there are, always accept the next one, and keep track of how many continuous
#timestamps there are
#get array of the distances between the two points
lc1=[]
slc=size(lc)
slc=slc[1]
for i in 1:slc-1
    push!(lc1,(lc[i+1]-lc[i]))
end
lc1=lc1/4
#find the regions of timesteps where there are continuous sequences of 4s, and how long they are
slc1=size(lc1)[1]
lc2=[]
for i in 1:slc1
    if lc1[i]==1
        push!(lc2,i)
    end
end
slc2=size(lc2)[1]
lc3=[]
for i in 1:slc2-1
    push!(lc3,(lc2[i+1]-lc2[i]))
end
lc3=convert(Array{Int16,1}, lc3)
lc4=lc3-ones(Int16, size(lc3))
