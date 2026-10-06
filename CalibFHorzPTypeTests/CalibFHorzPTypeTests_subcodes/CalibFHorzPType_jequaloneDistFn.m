function jequaloneDist=CalibFHorzPType_jequaloneDistFn(a_grid,z_grid,n_a,n_z)
% jequaloneDist as a function, used by CalibFHorzPTypeTests part 5: everyone starts with no assets at the middle shock (both types),
% which is exactly the matrix jequaloneDist of the bank's driver. The toolkit passes (a_grid,z_grid,n_a,n_z[,parameters]); a_grid and
% z_grid are not needed here.

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

end
