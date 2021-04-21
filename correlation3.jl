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



clusters7=[]
radius_1=[] # the inverse of distance between largest clusters

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
    clusters3=Array{Int64}(undef,2)
    clusters3[1]=clusters2[1]
    clusters3[2]=clusters2[2]
    clusters4=[]
    clusters5=[]
    for j in 1:2000
    if clusters[j]==clusters3[1]
        push!(clusters4,j)
    elseif clusters[j]==clusters3[2]
        push!(clusters5,j)
    end
    end
# we now want to calulate the average correlation between the largest clusters
# we ultimately want to graph these values as a function of time,N1N2,sqrt(N1N2),and 1/r.
# know time from the dataset, N1N2 functions from earlier in the code, and we can calculate 1/r by
# taking the center of mass of the clusters and finding the absolute difference

    #this is the calculation of the radius between large clusters
    c31=clusters3[1]
    c32=clusters3[2]
    # we have to define the center of mass for each cluster
    posx=g[SP[i]:EP[i],1]#x vector
    posy=g[SP[i]:EP[i],2]#y vector
    pos1x=[]
    pos1y=[]
    pos2x=[]
    pos2y=[]
    for p in 1:c31
    push!(pos1x,posx[clusters4[p]])
    end
    for p in 1:c31
    push!(pos1y,posy[clusters4[p]])
    end
    for p in 1:c32
    push!(pos2x,posx[clusters5[p]])
    end
    for p in 1:c32
    push!(pos2y,posy[clusters5[p]])
    end
    CM1=[mean(pos1x), mean(pos1y)]
    CM2=[mean(pos2x), mean(pos2y)]
    Distance=CM1-CM2
    rad=norm(Distance)
    rad1=1/rad
    push!(radius_1,rad1)
#correlation
# the correlation will be the average of the inner products between each pair of particles of different clusters
    velx=g[SP[i]:EP[i],4]
    vely=g[SP[i]:EP[i],5]

    clusters6=[] # we want to collect the correlations for each particle, but we want to ultimately reset the collection when we move to a new particle
    vel1=Array{Float64,2}(undef, c31, 2)
    vel2=Array{Float64,2}(undef, c32, 2)
    for k in 1:c31
    vel1[k,:]=[velx[clusters4[k]],vely[clusters4[k]]]
    end
    for l in 1:c32
    vel2[l,:]=[velx[clusters5[l]],vely[clusters5[l]]]
    end

    for m in 1:c31
        for o in 1:c32
            function f(m,o)
            dot(vel1[m,:],vel2[o,:])
            end
            push!(clusters6,f(m,o))
        end
    end
push!(clusters7,mean(clusters6))
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
