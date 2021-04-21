using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_Config.dat")
h=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_LargestCluster.dat")
# for a chosen timestep, plot the distribution of correlations between

#choose a timestep
#TS=300
n=2000
#create range of values to evaluate from timestep
g=f[38:end,1:9]
#starting point is at 300*number of particles
#end point is starting point plus number of particles



#the added term is due to there beeing a break between two timesteps, and
#this adding up over the space

# we now want to measure the velocity correlations between particles in different large
# clusters. we will measure the correlations between the biggest two clusters
# so we will need those values


e=h[37:end,1:8]
lp=j[:,1]
#clusters=g[SP:EP,7]
e1=e[:,1]
endtime=762
SP=Array{Int32}(undef, endtime)
EP=Array{Int32}(undef, endtime)
for x in 2:endtime

    SP[x]=(x-1)*(n+1)+2
    EP[x]=(x-1)*(n+1)+2+(n-1)
end




sizefcn1=[] # product of largest two clusters
sizefcn2=[] #product of largest two clusters, sqrt
sizefcn3=[]
sizefcn4=[]
for i in 2:endtime
    #function SP(x)
    #    (x-1)*(n+1)+2
    #end

    #function EP(x)
    #    (x-1)*(n+1)+2+(n-1)
    #end
    clusters=g[SP[i]:EP[i],7]
    clusters1=unique(clusters)
    clusters2=sort(clusters1, rev=true)
    clusters3=Array{Int32}(undef,2)
    clusters3[1]=clusters2[1]
    clusters3[2]=clusters2[2]

# we now want to calulate the average correlation between the largest clusters
# we ultimately want to graph these values as a function of time,N1N2,sqrt(N1N2),and 1/r.
# know time from the dataset, N1N2 functions from earlier in the code, and we can calculate 1/r by
# taking the center of mass of the clusters and finding the absolute difference

    #this is the calculation of the radius between large clusters
    c31=clusters3[1]
    c32=clusters3[2]
    product=c31*c32
    rat=c31/c32
    sqrat=sqrt(rat)
    sq=sqrt(product)
    push!(sizefcn1,product)
    push!(sizefcn2,sq)
    push!(sizefcn3,rat)
    push!(sizefcn4,sqrat)
end
