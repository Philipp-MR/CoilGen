function cut_circle=build_cut_circle(center_point,cut_witdh)
%build a rectangular cut shape

circular_resolution=10;

%build circular cut shapes
opening_circle=[sin(0:(2*pi)/(circular_resolution-1):2*pi); cos(0:(2*pi)/(circular_resolution-1):2*pi)];

% %create a circular opening cut
cut_circle=opening_circle.*repmat(cut_witdh/2,[2 1])+center_point;


end