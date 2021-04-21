# with this code we want to measure the radial distribution function
#between active matter clusters of defined sizes
#this is to measure whether there is a force between large particles.
using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
using StatsBase
using Combinatorics
using Distances
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_Config.dat")
g=f[38:end,1:9]
#the first series may need to have
#
sizeg=size(g)[1]
#start=660
#endtime=690
n=2000
endtime=Int64(sizeg/(n+1))
SP=Array{Int32}(undef, endtime)
EP=Array{Int32}(undef, endtime)

for x in 2:endtime

    SP[x]=(x-1)*(n+1)+2
    EP[x]=(x-1)*(n+1)+2+(n-1)
end
# we need to define the large clusters, and the different thresholds for measuring the radial distribution function
# for each timestep, define the large clusters, and define their positions. we can program in a method that does calculations above certain thresholds

# we will also need to define the clusters where particles belong in different halves/quadrants of the area
# due to boundary conditions, and we will need to

# defining the large clusters will be the first step, then we can separate them into smaller and smaller batches for different measurements.

state1=g[:,7]
xpos=g[:,1]
ypos=g[:,2]
 # these are the values that show which particles are in clusters of a given size. Index=particle, Value=size of cluster
statematrix=Matrix{Int64}(undef, 762, 2001) # the state matrix will be a matrix of n x m dimensions. n=timesteps  m=number of particles + 1(zeros row)
for i in 1:1524762      #for i in 1:(nm)
       j=Int64(ceil(i/2001))
       k=i-(2001*(j-1))
       statematrix[j,k]=state1[i]
end
statematrix=statematrix[:,2:end]
statematrix1=Matrix{Int64}(undef,761,2000)
for i in 1:761
    c=counts(statematrix[(i+1),:])
    siz=size(c)[1]
    for j in 1:siz
        statematrix1[i,j]=c[j]
    end
end


for i in 1:761
       for j in 1:2000
       statematrix1[i,j]=Int64(statematrix1[i,j]/j)
       end
end
N=100
statematrix2=statematrix1[:,N:2000] #statematrix2 gives the state matrix over the threshold
sm2size=size(statematrix2)
statematrix3=Matrix{Int64}(undef,1,(2000-N+1))
state3=[]
for i in 1:sm2size[1]
    if maximum(statematrix2[i,:])==1
        if size(findall(x->x==1,statematrix2[i,:]))[1]>1
        global statematrix3=vcat(statematrix3,transpose(statematrix2[i,:]))
        push!(state3,i)
        end
    end
end
statematrix3=statematrix3[2:end,:]
#now we have these timesteps and the statematrix corresponding to where there are only 1 of each large cluster size
#now we can begin calculating the radial distribution function
sm3size=size(state3)
resolution=.1
N=zeros(10000)
for i in 1:sm3size[1]


    #define the combinations of large particles that we have to measure
    sm3timestep=findall(x->x==1, statematrix3[i,:]) #find the sizes of the large clusters
    sm3timestep1=combinations(sm3timestep,2)
    sm3timestep1=collect(sm3timestep1)
    sm3timestepsize=size(sm3timestep)[1]
    sm3timestep1size=size(sm3timestep1)
     #if there's more than 1 large cluster # this has already been done
        for j in 1:sm3timestep1size[1]
            sizecluster1=sm3timestep1[j][1]#size of 1st cluster we are looking at
            sizecluster2=sm3timestep1[j][2]#size of 2nd cluster we are looking at
            #now we want to know which particles at this timestep that we are at, are in each of the clusters
            #we want to connect the index from the config file to a
            #for i in SP[state3[i]]:EP[state3[i]]
            indices1=findall(x->x==(sizecluster1+99),state1[SP[state3[i]+1]:EP[state3[i]+1]])
            indices2=findall(x->x==(sizecluster2+99),state1[SP[state3[i]+1]:EP[state3[i]+1]])
            distances=[] #empty array for the distances beteween particles, this will reset after each cluster pair is evaluated
            for k in 1:(sizecluster1+99)  # now we will calculate the distances between the particles
                for l in 1:(sizecluster2+99) #then we will determine the
                        pos1=[xpos[indices1[k]+SP[state3[i]+1]-1],ypos[indices1[k]+SP[state3[i]+1]-1]]
                        pos2=[xpos[indices2[l]+SP[state3[i]+1]-1],ypos[indices2[l]+SP[state3[i]+1]-1]]
                        push!(distances,euclidean(pos1,pos2))
                end
            end
            mdc=minimum(distances)
            mdc1=mdc/resolution
            mdc2=ceil(Int64,mdc)
            N[mdc2]=N[mdc2]+1
        end
end

N1=zeros(10000)
for i in 1:10000
    N1[i]=N[i]/(i*.1)
end
