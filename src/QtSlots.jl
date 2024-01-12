
using Blink
using NIfTI
using LinearAlgebra, StatsBase


plot3D(gx,gy,gz,Δf) = begin
    d = 0.0;        # Distance from (0,0,0)
    w = 0.0;        # Width

    n_slices = 1    # Number of slices

    downsample_factor = 6;


    # d1 = (d-w/2)*sqrt(gx^2+gy^2+gz^2);
    # d2 = (d+w/2)*sqrt(gx^2+gy^2+gz^2);

    traces = GenericTrace[];

    path = @__DIR__; path = path*"/datatypes/phantom/brain004/T1map.nii.gz";

    # Phantom
    ni = niread(path);
    data = ni.raw;

    # Clip outliers
    valorPercentil = percentile(data[:],99);
    data[data.>=valorPercentil].=valorPercentil;

    # Downsampling
    dsample_data = (data[1:downsample_factor:end,1:downsample_factor:end,1:downsample_factor:end]);

    # Normalization
    mini, maxi = extrema(dsample_data)
    normed = Float32.((dsample_data .- mini) ./ (maxi - mini));

    # # Phantom saturarion
    # sat = 0.05;  # Saturarion value (between 0 and 1). If there is so much difference between data, this value must be small, so we can see the all the range
    # V = dsample_data.>sat;
    # W = 1 .-V;
    # dsample_data .*= W;
    # dsample_data .+= (V.*sat);

	M, N, L = size(dsample_data)

    # FOV in mm
	FOVx = (M-1)*downsample_factor*1e-1
	FOVy = (N-1)*downsample_factor*1e-1
	FOVz = (L-1)*downsample_factor*1e-1

    x = -size(dsample_data)[1]/2 : (size(dsample_data)[1]/2) -1;
    y = -size(dsample_data)[2]/2 : (size(dsample_data)[2]/2) -1;
    z = -size(dsample_data)[3]/2 : (size(dsample_data)[3]/2) -1;

    ratio = FOVx/M

    X,Y,Z = mgrid(x,y,z);

    phantom = volume(x=X[:].*ratio,
                     y=Y[:].*ratio,
                     z=Z[:].*ratio,
                     value = normed[:],
                     opacity=0.1,
                     surface_count= 8,
                     isomin= 0.1,
                     isomax= 0.9,
                     # showscale = false,
                     colorscale = colors.RdBu_3
    )

    push!(traces,phantom);

    # Slices
    X,Y = mgrid(x,y);
    Z = 0*X;

    op = 0.3;

    # if n_slices > 1
    #     v = 0:n_slices-1;
    #     d = ((d-w/2) .+ v.*(w./(n_slices-1))).*sqrt(gx^2+gy^2+gz^2);
    # elseif n_slices == 1
    #     d = d .*sqrt(gx^2+gy^2+gz^2);
    # end

    d = 100 * Δf/γ

    for i in 1:n_slices
        if gz == 0
            if gy == 0
                trace = surface(x=(d[i]/gx).+Z,y=Y,z=X,
                                opacity=op,
                                showscale = false);
            else
                trace = surface(x=X,y=(d[i].-X*gx)/gy,z=Y,
                                opacity=op,
                                showscale = false);
            end
        else
            trace = surface(x=X,y=Y,z=(d[i].-X*gx.-Y*gy)/gz,
                            opacity=op,
                            showscale = false);
        end
        push!(traces,trace);
    end

    # Layout and plotting
    layout = Layout(
        scene_xaxis_range = [-FOVx/2,FOVx/2],
        scene_yaxis_range = [-FOVy/2,FOVy/2],
        scene_zaxis_range = [-FOVz/2,FOVz/2],
        scene_aspectratio = attr(x=.02*size(dsample_data)[1],
                                 y=.02*size(dsample_data)[2],
                                 z=.02*size(dsample_data)[3])
    )

    w = Window();
    p = plot(traces,layout);
    body!(w,p);

    # Necessary to see the plot from C++
    ww = Window();
    close(ww);
end




sim(mat,vec) = begin

    let
        #Window
        global w = nothing

        #Phantom
        ss = 5
#        phantom = KomaMRI.nii_brain_phantom3D(ss)
        phantom = KomaMRI.brain_phantom3D()

        #Scanner
        sys = Scanner()
        sys.B0 =        vec[1]       # Main magnetic field [T]
        sys.B1 =        vec[2]       # Max RF amplitude [T]
        sys.ADC_Δt =    vec[3]       # ADC sampling time
        sys.Gmax =      vec[4]       # Max Gradient [T/m]
        sys.Smax =      vec[5]       # Max Slew-Rate

        #Secuence
        global seq = Sequence()

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
                            n[3]*n[1]*(1-cos(θ))-n[2]*sin(θ)    n[3]*n[2]*(1-cos(θ))+n[1]*sin(θ)     cos(θ)+n[3]^2*(1-cos(θ)) ];


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

                global N_x = N
             end
         end

        seq.DEF = Dict("Nx"=>N_x,"Ny"=>N_x,"Nz"=>1)

        display(seq)

         #Simulation
 		global simParams = Dict{String,Any}()

        signal = simulate(phantom, seq, sys; simParams, w)
		signal

        ## Reconstruction
        # Nx, Ny = [N_x;N_x]
        # kdata = reshape(signal,(Nx,Ny))
        # kdata[:,2:2:Ny] = kdata[Nx:-1:1,2:2:Ny]
        # kdata = convert(Array{Complex{Float64},2},kdata)
		#
        # kdata_abs = abs.(kdata);
        # image = abs.(fftshift(ifft(kdata)))
		#
        # percentileValue = percentile(kdata_abs[:],99);
        # kdata_abs[kdata_abs.>=percentileValue].=percentileValue;


        # PLOT RECONSTRUCTED IMAGE FROM JULIA (NOT USED BECAUSE PLOTTING IS MADE IN C++)
        # hh, ww = 420,420

        # l = PlotlyJS.Layout(;title="Reconstruction",yaxis_title="y",
        #                     xaxis_title="x",height=hh,width=ww,
        #                     yaxis=attr(scaleanchor="x"),
        #                     modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest"
        #                     )

        # PlotlyJS.plot(PlotlyJS.heatmap(z=image,showscale=false,colorscale="Greys",transpose=true),l)


        # DISPLAY HISTOGRAM OF RECONSTRUCTED IMAGE --------------------------------------------
    #    l = PlotlyJS.Layout(barmode="overlay")
    #    plt = PlotlyJS.plot(PlotlyJS.histogram(x=vec(abs_ifft)),l)

    #    w = Window();
    #    body!(w,plt);

    #    # Necessary to see the plot from C++
    #    ww = Window();
    #    close(ww);
        # -------------------------------------------------------------------------------------

        # [image kdata_abs]
    end
end

# This function returns to c++ the final sequence in a matrix-form, so that we are able to plot it on qml
displaySeq() = begin
    aux = zeros(6,length(seq))

    for i=1:length(seq)
        aux[1,i] = seq[i].DUR[1]                                                          #duration
        aux[2,i] = abs(seq[i].RF[1].A)>0.0 ? (seq[i].GR[1].A/100) : (seq[i].GR[1].A)      #gx
        aux[3,i] = abs(seq[i].RF[1].A)>0.0 ? (seq[i].GR[2].A/100) : (seq[i].GR[2].A)      #gy
        aux[4,i] = abs(seq[i].RF[1].A)>0.0 ? (seq[i].GR[3].A/100) : (seq[i].GR[3].A)      #gz
        aux[5,i] = seq[i].RF[1].A                                                         #RF amplitude
        aux[6,i] = seq[i].ADC[1].N>0 ? 1 : 0                                              #Readout on/off
    end

    aux
end
