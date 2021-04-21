using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_Config.dat")
h=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_LargestCluster.dat")
# for a chosen timestep, plot the distribution of correlations between

#choose a timestep
TS=300
n=2000
#create range of values to evaluate from timestep
g=f[38:end,1:9]
#starting point is at 300*number of particles
#end point is starting point plus number of particles
SP=(TS-1)*(n+1)+2
#the added term is due to there beeing a break between two timesteps, and
#this adding up over the space
EP=SP+(n-1)
# we now want to measure the velocity correlations between particles in different large
# clusters. we will measure the correlations between the biggest two clusters
# so we will need those values
j=h[37:end,1:8]
lp=j[:,1]
clusters=g[SP:EP,7]
clusters1=unique(clusters)
clusters2=sort(clusters1, rev=true)
clusters3=clusters2[1:2]
clusters4=[]
clusters5=[]
for i in 1:2000
    if clusters[i]==clusters3[1]
        push!(clusters4,i)
    elseif clusters[i]==clusters3[2]
        push!(clusters5,i)
    end
end
# we now want to calulate the average correlation between the largest clusters
# we ultimately want to graph these values as a function of time,N1N2,sqrt(N1N2),and 1/r.
# know time from the dataset, N1N2 functions from earlier in the code, and we can calculate 1/r by
# taking the center of mass of the clusters and finding the absolute difference

#correlation
# the correlation will be the average of the inner products between each pair of particles of different clusters
velx=g[SP:EP,4]
vely=g[SP:EP,5]
c31=clusters3[1]
c32=clusters3[2]
clusters6=Array{Float64}(undef, (c31*c32) ) # we want to collect the correlations for each particle, but we want to ultimately reset the collection when we move to a new particle
vel1=Array{Float64,2}(undef, c31, 2)
vel2=Array{Float64,2}(undef, c32, 2)
for i in 1:c31
    vel1[i,:]=[velx[clusters4[i]],vely[clusters4[i]]]
end
for j in 1:c32
    vel2[j,:]=[velx[clusters5[j]],vely[clusters5[j]]]
end

for i in 1:c31
    for j in 1:c32
        function f(i,j)
            dot(vel1[i,:],vel2[j,:])
        end
        push!(clusters6,f(i,j))
    end
end



# we now have our particles we want to look at
#we need to establish which are which, so then we can evaluate correlation functions
#clusters6=[]
#clusters7=[]
#for i in 1:2000
#    if clusters[i]==clusters3[1]
#        push!(clusters6,1)
#    elseif clusters[i]==clusters3[2]
#        push!(clusters7,2)
#    end
#end
# we now want to find the distance between the particles, as well as the inner
#product of the velocities
#we want the radius of the
