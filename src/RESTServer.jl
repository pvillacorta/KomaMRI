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

using StatsBase, FFTW
using StructTypes

@everywhere begin
   using SharedArrays
   using KomaMRI
   using LinearAlgebra
   using JSON3
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
@get "/greet" function(req::HTTP.Request)
   return "Hello world!"
end

## SIMULATION
@post "/simulate" function(req::HTTP.Request)
   global aux = JSON3.read(req.body)
   N_blocks = length(aux.mat[1])     # Number of columns (blocks)

   mat::Matrix{Float64} = zeros(10,N_blocks)
   for i in range(1,10,step=1)
      mat[i,:] = aux.mat[i]
   end

   vec::Vector{Float64} = aux.vec
   phantom::String = aux.phantom

   global statusFile = tempname()
   touch(statusFile)
   print(statusFile, '\n')

   # Simulation  (asynchronous. It should not block the HTTP 202 Response)
   global result = @spawnat 2 sim(mat,vec,phantom,statusFile)          # Process 2 executes simulation
   # global result = remotecall(sim, 2, mat, vec, statusFile)      # Equivalent to expression above

   # while 1==1
   #    io = open(statusFile,"r")
   #    if (!eof(io))
   #       global simProgress = read(io,Int32)
   #       print("leido\n")
   #    end
   #    close(io)
   #    print("Progreso: ", simProgress, '\n')
   #    sleep(0.2)
   # end

   # Actualizar tabla de correspondencias simulación-proceso

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
      global simProgress = read(io,Int32)
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


## PLOT SEQUENCE
# @post "/plot" function(req::HTTP.Request)
#    data = json(req)
#    sys = Scanner()
#    json_to_seq(data,sys)
#    return HTTP.Response(200,body="Exito")
# end
# ---------------------------------------------------------------------------

serve(host="0.0.0.0",port=8085)
