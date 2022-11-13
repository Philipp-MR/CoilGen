function input=temp_evaluation(input,target_field)
%Evaluate wether precalculated values can be used from previous
%calculations

preoptimization_input_hash=generate_DataHash({input.coil_mesh_file input.iteration_num_mesh_refinement input.surface_is_cylinder_flag,target_field});
optimized_input_hash=generate_DataHash({input.sf_opt_method input.tikonov_reg_factor input.fmincon_parameter});

input.temp_evalution.use_preoptimization_temp=false;
input.temp_evalution.use_optimized_temp=false;

%initialize values if not exisiting
if ~isfield(input.temp,'preoptimization_hash')
input.temp.preoptimization_hash='none';
end
if ~isfield(input.temp,'optimized_hash')
input.temp.optimized_hash='none';
end

if strcmp(preoptimization_input_hash,input.temp.preoptimization_hash)
input.temp_evalution.use_preoptimization_temp=true;
end

if strcmp(optimized_input_hash,input.temp.optimized_hash)
input.temp_evalution.use_optimized_temp=true;
end

%Assign the new hash to temp
input.preoptimization_hash=preoptimization_input_hash;
input.optimized_hash=optimized_input_hash;

function H = generate_DataHash(Data)
Engine = java.security.MessageDigest.getInstance('MD5');
H = CoreHash(Data, Engine);
H = sprintf('%.2x', H);   % To hex string
function H = CoreHash(Data, Engine)
% Consider the type of empty arrays:
S = [class(Data), sprintf('%d ', size(Data))];
Engine.update(typecast(uint16(S(:)), 'uint8'));
H = double(typecast(Engine.digest, 'uint8'));
if isa(Data, 'struct')
   n = numel(Data);
   if n == 1  % Scalar struct:
      F = sort(fieldnames(Data));  % ignore order of fields
      for iField = 1:length(F)
         H = bitxor(H, CoreHash(Data.(F{iField}), Engine));
      end
   else  % Struct array:
      for iS = 1:n
         H = bitxor(H, CoreHash(Data(iS), Engine));
      end
   end
elseif isempty(Data)
   % No further actions needed
elseif isnumeric(Data)
   Engine.update(typecast(Data(:), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif ischar(Data)  % Silly TYPECAST cannot handle CHAR
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif iscell(Data)
   for iS = 1:numel(Data)
      H = bitxor(H, CoreHash(Data{iS}, Engine));
   end
elseif islogical(Data)
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif isa(Data, 'function_handle')
   H = bitxor(H, CoreHash(functions(Data), Engine));
else
   warning(['Type of variable not considered: ', class(Data)]);
end
end
end

end