#here we want to look at the correlation functions of
#clusters that collide with each other.
using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
cd("C:/Users/AB/Documents/Computation/ClusterCode")
f=readdlm("Series_02_Angle_45.0_NP_2000_NV_4_Rho_0.200_Config.dat")
g=f[38:end,1:9]

start=660
endtime=690
SP=Array{Int32}(undef, endtime)
EP=Array{Int32}(undef, endtime)
n=2000
for x in 2:endtime

    SP[x]=(x-1)*(n+1)+2
    EP[x]=(x-1)*(n+1)+2+(n-1)
end
Length=200

C=zeros(10000)
N=zeros(10000)
Cavg=zeros(10000)
delta=.1
#C1=zeros((endtime-start+1))
#CRad=zeros(100000)
for i in start:endtime

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
        pos1=Array{Float64,2}(undef, c31, 2)
        pos2=Array{Float64,2}(undef, c32, 2)
        for k in 1:c31
        pos1[k,:]=[posx[clusters4[k]],posy[clusters4[k]]]
        end
        for l in 1:c32
        pos2[l,:]=[posx[clusters5[l]],posy[clusters5[l]]]
        end

    #correlation
    # the correlation will be the average of the inner products between each pair of particles of different clusters
        velx=g[SP[i]:EP[i],4]
        vely=g[SP[i]:EP[i],5]

        #clusters6=[] # we want to collect the correlations for each particle, but we want to ultimately reset the collection when we move to a new particle
        vel1=Array{Float64,2}(undef, c31, 2)
        vel2=Array{Float64,2}(undef, c32, 2)
        for k in 1:c31
        vel1[k,:]=[velx[clusters4[k]],vely[clusters4[k]]]
        end
        for l in 1:c32
        vel2[l,:]=[velx[clusters5[l]],vely[clusters5[l]]]
        end
rlist=[]
        for m in 1:c31
            for o in 1:c32
                #we want to calculate the distance between two particles
                #then we will determine which of our bins we adjust
                rminmat=zeros(3,3)
                for p in [-1,0,1]
                    for q in [-1,0,1]
                        function g(p,q)

                            (sqrt((((pos2[o,1]-pos1[m,1]+(p*Length))^(2))+(((pos2[o,2]-pos1[m,2]+(q*Length))^(2))))))
                        end
                        rminmat[(p+2),(q+2)]=g(p,q)
                    end
                end
                rmin=minimum(rminmat)
                vals=argmin(rminmat)
                p1=vals[1]
                q1=vals[2]
                r12=[(pos2[o,1]-pos1[m,1]+((p1-2)*Length)),(pos2[o,1]-pos1[m,1]+((q1-2)*Length))]
                r=norm(r12)

                push!(rlist,r)
                #function f(m,o)
                #    dot(r12,(vel1[m,:]-vel2[o,:]))/r
                #end
                #C1[(i-start+1)]=C1[(i-start+1)]+f(m,o)


            end
        end

                #C1[(i-start+1)]=C1[(i-start+1)]/(c31*c32)

                #rfloor=r/delta
                #rfloor1=ceil(Int,rfloor)
                #N[rfloor1]=N[rfloor1]+1
                #radtimei=minimum(rlist)
                #radtimei1=Int64(floor(radtimei/delta))
                #CRad[radtimei1]=C1[(i-start+1)]

                function f(m,o)
                    dot(r12,(vel1[m,:]-vel2[o,:]))/r
                end


                C[rfloor1]=C[rfloor1]+f(m,o)
                #function f(m,o)
                #dot(vel1[m,:],vel2[o,:])
                #end
            #end
        end

   for i in 1:10000
        Cavg[i]=C[i]/N[i]
    end
