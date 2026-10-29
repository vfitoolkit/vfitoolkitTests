function jequaloneDist=CalibSGE_jequaloneDistFn(a_grid,z_grid,n_a,n_z)
% jequaloneDist as a function, used by CalibStationaryGeneralEqmTests part 5: everyone starts with zero assets and z drawn from the
% stationary distribution of the bank's z process, which is exactly the matrix jequaloneDist of the bank's driver. The toolkit passes
% (a_grid,z_grid,n_a,n_z[,parameters]); a_grid and z_grid are not needed here.

[~,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
pi_z_stat=pi_z^1000; pi_z_stat=pi_z_stat(1,:)';
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,:)=shiftdim(pi_z_stat,-1);

end
