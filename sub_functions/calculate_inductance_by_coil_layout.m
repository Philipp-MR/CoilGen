function coil_parts =calculate_inductance_by_coil_layout(coil_parts,input)
%Calculate the inducate by means of integration of the vector potential
%along the coil track


conductor_width=input.conductor_cross_section_width;
conductor_height=input.conductor_cross_section_height;
skip_inductance_calculation=input.skip_inductance_calculation;


down_sample_factor=10; % reduce the number of track points since fasthenry has an upper limit..

coil_parts(numel(coil_parts)).coil_resistance=[];
coil_parts(numel(coil_parts)).coil_inductance=[];
coil_parts(numel(coil_parts)).coil_length=[];
coil_parts(numel(coil_parts)).coil_cross_section=[];



if ~skip_inductance_calculation


sim_freq=1; %in Hz
material_conductivity=5.8e7; % in 1/(mm*Ohm)

try
[Result, ~] = system('WHERE /F /R "c:\Program Files (x86)" FastHenry2.exe');
catch
Result=1;
end

if Result==0
    
for part_ind=1:numel(coil_parts)
[coil_parts(part_ind).coil_resistance,coil_parts(part_ind).coil_inductance,coil_parts(part_ind).coil_cross_section]=create_fast_henry_file(coil_parts(part_ind).wire_path,conductor_width,conductor_height,sim_freq,material_conductivity,down_sample_factor);
end

else
    fprintf(' FastHenry2 is not installed in the Folder- "Program Files (x86)" -", magnetic inductance will not be calculated.. \n');
    
    for part_ind=1:numel(coil_parts)
    coil_parts(part_ind).coil_resistance=0;
    coil_parts(part_ind).coil_inductance=0;
    coil_parts(part_ind).coil_length=0;
    coil_parts(part_ind).coil_cross_section=0;
    end
end

else
    for part_ind=1:numel(coil_parts)
    coil_parts(part_ind).coil_resistance=0;
    coil_parts(part_ind).coil_inductance=0;
    coil_parts(part_ind).coil_length=0;
    coil_parts(part_ind).coil_cross_section=0;
    end
    
end


%Calculate the length of the coil
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).coil_length=sum(vecnorm(coil_parts(part_ind).wire_path.v(:,2:end)-coil_parts(part_ind).wire_path.v(:,1:end-1))); % in m
end

function [coil_resistance,coil_inductance,coil_cross_section]=create_fast_henry_file(wire_path,conductor_width,conductor_height,sim_freq,material_conductivity,down_sample_factor)
%Create a FASTHENRY2 input file, run it with fasthenry2 and read out the result


    %INPUTS:
    %wire_path: sequence of 3dim points coords in m, Matrix: 3x(Number of points)
    %conductor_width: width of conductor in m
    %conductor_height: height of conductor in m
    %sim_freq: in hz
    %material_conductivity: conductivity in 1/(mm*ohm)

%     conductor_width=0.0001;
%     conductor_height=0.0001;
    

wire_path_downsampled.uv=wire_path.uv(:,1:down_sample_factor:end);
wire_path_downsampled.v=wire_path.v(:,1:down_sample_factor:end);

%Create the ".inp" input file
        fast_henry_file_name='coil_track_FH2_input.inp';
        %delete(fast_henry_file_name);
        fid = fopen(fast_henry_file_name, 'w+' );
        fprintf(fid,'\n\n\n.Units m\n'); %specify the unit system
        fprintf(fid,".Default sigma="+num2str(material_conductivity)+" \n\n");
        %define the points (vertices) in the FASTHENRY input file
        for point_ind = 1:size(wire_path_downsampled.v,2)
        str_to_include="\n"+"N"+num2str(point_ind)+" "+"x="+num2str(wire_path_downsampled.v(1,point_ind))+" "+"y="+num2str(wire_path_downsampled.v(2,point_ind))+" "+"z="+num2str(wire_path_downsampled.v(3,point_ind));
        fprintf(fid,str_to_include);
        end
        fprintf(fid,'\n\n\n');
        %define the cross section of the wire segment
        for seg_ind = 2:size(wire_path_downsampled.v,2)
        str_to_include="\n"+"E"+num2str(seg_ind-1)+" "+ "N"+num2str(seg_ind-1) +" "+ "N"+num2str(seg_ind)+" " +"w="+num2str(conductor_width)+" " +"h="+num2str(conductor_height);
        fprintf(fid,str_to_include);
        end
        fprintf(fid,"\n\n"+".external N1 N"+num2str(size(wire_path_downsampled.v,2)));
        fprintf(fid,"\n\n"+".freq fmin="+" "+num2str(sim_freq) +" "+" fmax="+num2str(sim_freq)+" "+ "ndec=1");
        fprintf(fid,"\n\n"+".end");
        fclose(fid);
    
%Create the windows ".vbs" script file to run the ".inp" automatically
        script_file_name='run_FH2.vbs';
        fid2 = fopen(script_file_name, 'w+' );
        fprintf(fid2,'Set FastHenry2 = CreateObject("FastHenry2.Document")\n');
        fprintf(fid2,'%s',string(strcat('couldRun = FastHenry2.Run("',cd,'\',fast_henry_file_name,'")')));
        fprintf(fid2,'\nDo While FastHenry2.IsRunning = True');
        fprintf(fid2,'\nWscript.Sleep 500');
        fprintf(fid2,'\nLoop');
        fprintf(fid2,'\ninductance = FastHenry2.GetInductance()');
        fprintf(fid2,'\nFastHenry2.Quit');
        fprintf(fid2,'\nSet FastHenry2 = Nothing');
        fclose(fid2);

        
%Run the script
        
system(char(string(cd)+"\"+script_file_name));


%Read the results

fid3=fopen('Zc.mat');
out=textscan(fid3,'%s','headerlines',2);
fclose(fid3);

real_Z=str2num(out{1}{1});
im_Z=imag(str2num(out{1}{2}));

coil_cross_section=conductor_width*conductor_height; % in m²

coil_inductance=im_Z/(sim_freq*2*pi); %in Henry
coil_resistance=real_Z; % in Ohm

%Remove the created script files again
delete Zc.mat
delete run_FH2.vbs
delete coil_track_FH2_input.inp

end




end