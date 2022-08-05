function track_out= smooth_track_by_folding(track_in,smoothing_length)
%@Philipp Amrein, 2022 , Uniklinik Freiburg


track_out=track_in;

if smoothing_length>0

extended_track=[repmat(track_in(:,1),[1 smoothing_length]) track_in(:,2:end-1) repmat(track_in(:,end),[1 smoothing_length])];
for shift_ind=[(-1)*(smoothing_length-1):-1 1:(smoothing_length-1)]
add_track=circshift(extended_track,shift_ind,2);
add_track=add_track(:,smoothing_length:(end-smoothing_length+1));
track_out=track_out+add_track;
end

track_out=track_out./(2*smoothing_length-1);

end

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

