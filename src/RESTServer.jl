using Distributed

if nprocs() <= 1
   addprocs(1)
end

@everywhere begin
   cd("/home/export/personal/pvilayl/KomaMRI.jl")
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
end

using Oxygen
using HTTP
using JSON3

using StatsBase, FFTW
using StructTypes

@everywhere begin
   using SharedArrays
   using KomaMRI
   using LinearAlgebra
end


staticfiles("content", "static")  # Give access to files in content directory
                                  # to access to them from the client, we should make
                                  # a request to http://serverdir:port/static/filename

# ------------------------------- STRUCTS ------------------------------------
global simulationId = 1
global simProgress = -1
global statusFile = ""

mutable struct Image
   # id::Int
   data::Array{Float32,2}
end

# ------------------------------ FUNCTIONS ------------------------------------
@everywhere include("ServerFunctions.jl")


# ---------------------------- API METHODS ---------------------------------
@get "/foo" function (req::HTTP.Request)
   global statusFile = tempname()
   touch(statusFile)

   mat  = [1        2      5;           # cod
           5.87e-4  0.01   0;           # dur
           0        0      0;           # gx
           0        0      0;           # gy
           1        0      0;           # gz
           10e-6    0      0;           # b1x
           0        0      0;           # b1y
           0        0      0;           # Δf
           0        0      0.4;         # fov
           0        0      201]         # n

   vec  =  [1.5;          #B0
            10e-6;        #B1
            2e-6;         #Delta_t
            60e-3;        #Gmax
            500]          #Smax

   image =  @spawnat 2 sim(mat, vec, statusFile)
   # image =  sim(mat, vec, statusFile)
end

@post "/simulate" function(req::HTTP.Request)
   global aux = JSON3.read(req.body)
   N_blocks = length(aux.mat[1])     # Number of columns (blocks)

   mat::Matrix{Float64} = zeros(10,N_blocks)
   for i in range(1,10,step=1)
      mat[i,:] = aux.mat[i]
   end

   vec::Vector{Float64} = aux.vec

   global statusFile = tempname()
   touch(statusFile)

   # Simulation  (asynchronous. It should not block the HTTP 202 Response)
   # global result = @spawnat 2 sim(mat,vec,statusFile)          # Process 2 executes simulation
   global result = remotecall(sim, 2, mat, vec, statusFile)      # Equivalent to expression above

   # while 1==1
   #    io = open(statusFile,"r")
   #    if (!eof(io))
   #       global simProgress = read(io,Int64)
   #    end
   #    close(io)
   #    print("Progreso: ", simProgress, '\n')
   #    sleep(0.2)
   # end


   headers = ["Location" => string("/simulate/",simulationId)]
   global simulationId += 1
   # 202: Partial Content
   return HTTP.Response(202,headers)
end

"""
                  [ -1,      if the simulation has not started yet
    simProgress = [ (0,100), if the simulation is running
                  [ 100,     if the simulation has finished, reconstructing
                  [ 101,     if the reconstruction has finished
"""


"""If the simulation has finished, it returns its result. If not, it returns 303 with location = /simulate/{simulationId}/status"""

@get "/simulate/{simulationId}" function(req::HTTP.Request, simulationId)
   io = open(statusFile,"r")
   if (!eof(io))
      global simProgress = read(io,Int64)
   end
   close(io)

   if simProgress < 101      # Simulation not started or in progress
      headers = ["Location" => string("/simulate/",simulationId,"/status")]
      return HTTP.Response(303,headers)
   elseif simProgress == 101  # Simulation finished
      global simProgress = -1
      im = Image(fetch(result))
      return HTTP.Response(200,body=JSON3.write(im))
   end
end


@get "/simulate/{simulationId}/status" function(req::HTTP.Request, simulationId)
   return HTTP.Response(200,body=JSON3.write(simProgress))
end


# ---------------------------------------------------------------------------

serve(host="0.0.0.0",port=8085)
