function jequaloneDist=CalibFHorz_jequaloneDistFn(a_grid,z_grid,n_a,n_z)
% jequaloneDist as a function, used by CalibFHorzTests part 5: everyone starts with no assets at the middle shock, which is exactly
% the matrix jequaloneDist of the bank's driver. The toolkit passes (a_grid,z_grid,n_a,n_z[,parameters]); a_grid and z_grid are
% not needed here (z_grid arrives as the user's grid or the internal joint grid depending on simoptions.jequaloneDist_usergrids).

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

end
