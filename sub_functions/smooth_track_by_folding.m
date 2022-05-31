function track_out= smooth_track_by_folding(track_in,smoothing_length)
%@Philipp Amrein, 2022 , Uniklinik Freiburg

track_out=zeros(size(track_in));
%track_out.v=zeros(size(track_in.v));
for shift_ind=-smoothing_length:1:smoothing_length
track_out=track_out+track_in(:,circshift(1:size(track_in,2),shift_ind));
%track_out.v=track_out.v+track_in.v(:,circshift(1:size(track_in.v,2),shift_ind));
end

track_out=track_out./(2*smoothing_length+1);
%track_out.v=track_out.v./(2*smoothing_length+1);



end

