function [ vec ] = convertToVec( planar_pd, z_ang )


    vec = [cos(z_ang).*cos(planar_pd), cos(z_ang).*sin(planar_pd), sin(z_ang)];

end

