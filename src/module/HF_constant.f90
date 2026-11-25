module HF_parameter
    integer         ,parameter    ::   maxSCF = 100
    double precision,parameter    ::   thresh = 1d-9
    integer         ,parameter    :: max_diis = 8
    double precision, allocatable, save :: hf_P_last(:,:)   ! [Nov 2026] converged density (keyword fs_correction)
end module HF_parameter