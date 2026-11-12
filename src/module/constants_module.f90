module constants_module
    integer         ,parameter     :: max_atom = 1500
    double precision,parameter     :: pi       = acos(-1.d0)
    double precision,parameter     :: Ang_par  = 1.8897261249935897D00
    integer, parameter             :: dpc      = kind((1.0d0,1.0d0))
    integer, parameter             :: dp       = kind(1.0d0)
end module constants_module