function track_out= smooth_track_by_folding_within_3d(track_in,smoothing_length)
%@Philipp Amrein, 2022 , Uniklinik Freiburg

%IDEA: If the mesh has a strong curvature, the smooting within the 2D
%domian (merging neighbouring 2D points) will lead to directional errors

%if the smooting is done for 3D points, the might not be on the surface any
%more

%therefore are correction must be done afterwards where the 3D points are
%properly put on the surface , only then are the corresponding 2D point
%well defined



track_out=track_in;

% if smoothing_length>0
% 
% extended_track=[repmat(track_in(:,1),[1 smoothing_length]) track_in(:,2:end-1) repmat(track_in(:,end),[1 smoothing_length])];
% for shift_ind=[(-1)*(smoothing_length-1):-1 1:(smoothing_length-1)]
% add_track=circshift(extended_track,shift_ind,2);
% add_track=add_track(:,smoothing_length:(end-smoothing_length+1));
% track_out=track_out+add_track;
% end
% 
% track_out=track_out./(2*smoothing_length-1);
% 
% end

%%%Old version
% track_out=zeros(size(track_in));
% %track_out.v=zeros(size(track_in.v));
% for shift_ind=-smoothing_length:1:smoothing_length
% track_out=track_out+track_in(:,circshift(1:size(track_in,2),shift_ind));
% %track_out.v=track_out.v+track_in.v(:,circshift(1:size(track_in.v,2),shift_ind));
% end
% 
% track_out=track_out./(2*smoothing_length+1);
% %track_out.v=track_out.v./(2*smoothing_length+1);



end

