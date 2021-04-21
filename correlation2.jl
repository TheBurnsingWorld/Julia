using DelimitedFiles
using Plots
using LinearAlgebra
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
#SP=(TS-1)*(n+1)+2
#the added term is due to there beeing a break between two timesteps, and
#this adding up over the space
#EP=SP+(n-1)
# we now want to measure the velocity correlations between particles in different large
# clusters. we will measure the correlations between the biggest two clusters
# so we will need those values
j=h[37:end,1:8]
j1=j[:,1]
endtime=j1[762]
clusters7=Array{Float64}(undef,endtime)

for i in 1:endtime
           TS=i
           SP=(TS-1)*(n+1)+2
           EP=SP+(n-1)
           clusters=g[SP:EP,7]
           clusters1=unique(clusters)
           clusters2=sort(clusters1, rev=true)
           clusters3=clusters2[1:2]
           clusters4=[]
           clusters5=[]
           for j in 1:2000
               if clusters[j]==clusters3[1]
                   push!(clusters4,j)
               elseif clusters[j]==clusters3[2]
                   push!(clusters5,j)
               end
           end
           velx=g[SP:EP,4]
           vely=g[SP:EP,5]
           c31=clusters3[1]
           c32=clusters3[2]
           clusters6=Array{Float64}(undef, (c31*c32) ) # we want to collect the correlations for each particle, but we want to ultimately reset the collection when we move to a new particle
               vel1=Array{Float64,2}(undef, c31, 2)
           vel2=Array{Float64,2}(undef, c32, 2)
           for k in 1:c31
               vel1[k,:]=[velx[clusters4[k]],vely[clusters4[k]]]
                   end
           for l in 1:c32
               vel2[l,:]=[velx[clusters5[l]],vely[clusters5[l]]]
                   end

           for m in 1:c31
               for n in 1:c32
                   function f(i,j)
                       dot(vel1[m,:],vel2[n,:])
                   end
                   push!(clusters6,f(m,n))
               end
           end

           clusters7[i]=mean(clusters6)
end
