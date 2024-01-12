"Convert a 1D vector with system paramaters into a KomaMRICore.Scanner object"
vec_to_scanner(vec) = begin
   sys = Scanner()
   sys.B0 =        vec[1]       # Main magnetic field [T]
   sys.B1 =        vec[2]       # Max RF amplitude [T]
   sys.ADC_Δt =    vec[3]       # ADC sampling time
   sys.Gmax =      vec[4]       # Max Gradient [T/m]
   sys.Smax =      vec[5]       # Max Slew-Rate

   sys
end

"Convert a 2D matrix containing sequence information into a KomaMRICore.Sequence object"
mat_to_seq(mat,sys::Scanner) = begin
   seq = Sequence()

   for i=1:size(mat)[2]

      if mat[1,i] == 1 # Excitation
         B1 = mat[6,i] + mat[7,i]im;
         duration = mat[2,i]
         Δf = mat[8,i];
         EX = PulseDesigner.RF_hard(B1, duration, sys; G = [mat[3,i] mat[4,i] mat[5,i]], Δf)
         seq += EX

         G = EX.GR.A; G = [G[1];G[2];G[3]];
         REF = [0;0;1];

         # We need to create a rotation matrix which transfomrs vector [0 0 1] into vector G
         # To do this, we can use axis-angle representation, and then calculate rotation matrix with that
         # https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_rotation_matrix_to_axis%E2%80%93angle

         # Cross product:
         global cross_prod = LinearAlgebra.cross(REF,G);
         # Rotation axis (n = axb) Normalized cross product:
         n = normalize(cross_prod);
         # Rotation angle:
         θ = asin(norm(cross_prod)/((norm(REF))*(norm(G))));
         # Rotation matrix:
         global R = [cos(θ)+n[1]^2*(1-cos(θ))            n[1]*n[2]*(1-cos(θ))-n[3]*sin(θ)     n[1]*n[3]*(1-cos(θ))+n[2]*sin(θ);
                     n[2]*n[1]*(1-cos(θ))+n[3]*sin(θ)    cos(θ)+n[2]^2*(1-cos(θ))             n[2]*n[3]*(1-cos(θ))-n[1]*sin(θ);
                     n[3]*n[1]*(1-cos(θ))-n[2]*sin(θ)    n[3]*n[2]*(1-cos(θ))+n[1]*sin(θ)     cos(θ)+n[3]^2*(1-cos(θ))        ];


      elseif mat[1,i] == 2 # Delay
         DELAY = Delay(mat[2,i])
         seq += DELAY


      elseif ((mat[1,i] == 3) || (mat[1,i] == 4)) # Dephase or Readout
         AUX = Sequence()

         ζ = abs(sum([mat[3,i],mat[4,i],mat[5,i]])) / sys.Smax
         ϵ1 = mat[2,i]/(mat[2,i]+ζ)

         AUX.GR[1] = Grad(mat[3,i],mat[2,i],ζ)
         AUX.GR[2] = ϵ1*Grad(mat[4,i],mat[2,i],ζ)
         AUX.GR[3] = Grad(mat[5,i],mat[2,i],ζ)

         AUX.DUR = convert(Vector{Float64}, AUX.DUR)
         AUX.DUR[1] = mat[2,i] + 2*ζ

         AUX.ADC[1].N = (mat[1,i] == 3) ? 0 : trunc(Int,mat[10,i])     # No samples during Dephase interval
         AUX.ADC[1].T = mat[2,i]                                       # The duration must be explicitly stated
         AUX.ADC[1].delay = ζ

         AUX = (norm(cross_prod)>0) ? R*AUX : AUX

         if(mat[1,i]==4)
             global N_x = trunc(Int,mat[10,i])
         end

         seq += AUX


      elseif mat[1,i] == 5 # EPI
          FOV = mat[9,i]
          N = trunc(Int,mat[10,i])

          EPI = PulseDesigner.EPI(FOV, N, sys)
          EPI = (norm(cross_prod)>0) ? R*EPI : EPI
          seq += EPI

          N_x = N
      end

   end

   seq.DEF = Dict("Nx"=>N_x,"Ny"=>N_x,"Nz"=>1)

   seq
end

"Obtain the reconstructed image from raw_signal (obtained from simulation)"
recon(raw_signal) = begin
   recParams = Dict{Symbol,Any}(:reco=>"direct")

   acqData = AcquisitionData(raw_signal)
   acqData.traj[1].circular = false #Removing circular window
   acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqData.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
   Nx, Ny = raw_signal.params["reconSize"][1:2]
   recParams[:reconSize] = (Nx, Ny)
   recParams[:densityWeighting] = true

   aux = @timed reconstruction(acqData, recParams)
   image  = reshape(aux.value.data,Nx,Ny,:)
   kspace = fftc(reshape(aux.value.data,Nx,Ny,:))

   # Conversion to uint8
   image_aux = abs.(image[:,:,1])
   uint8_image = round.(Int,255*(image_aux./maximum(image_aux)))

   uint8_image
end

"Obtain raw RM signal. Input arguments are a 2D matrix (sequence) and a 1D vector (system parameters)"
sim(mat,vec,path) = begin
   print("Prueba")

   # Phantom
   phant = KomaMRI.brain_phantom2D()

   # Scanner
   sys = vec_to_scanner(vec)

   # Sequence
   seq = mat_to_seq(mat,sys)

   # Simulation parameters
   simParams = Dict{String,Any}()

   # Simulation
   raw_signal = simulate(phant, seq, sys, path; simParams, w=nothing)

   # Reconstruction
   image = recon(raw_signal)

   if path !== nothing
         io = open(path,"w") # "w" mode overwrites last status value, even if it has not been read yet
         write(io,trunc(Int,101))
         close(io)
   end

   image
end
