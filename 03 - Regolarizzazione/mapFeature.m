function out = mapFeature(phi1, phi2, degree)
% MAPFEATURE Feature mapping function to polynomial features
%
%   MAPFEATURE(phi1, phi2, degree) maps the two input features
%   to polynomial features of degree "degree"
%
%   Returns a new feature array with more features, comprising of 
%   phi1, phi2, phi1.^2, phi2.^2, phi1*phi2, phi1*phi2.^2, etc..
%
%   Inputs phi1, phi2 must be the same size
%

out = ones(size(phi1(:,1)));
for i = 1:degree
    for j = 0:i
        out(:, end+1) = (phi1.^(i-j)).*(phi2.^j);
    end
end

end