import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

-- local gamma = 1.4
local gamma = 5.0/3.0
local gamma_m1 = gamma - 1.0
local onebygm1 = 1.0 / gamma_m1

__demand(__inline)
task get_pressure( rho : double, rhoe : double )
  return gamma_m1 * rhoe
end

__demand(__inline)
task get_internal_energy( rho : double, p : double)
  return onebygm1 * p
end

__demand(__inline)
task get_sos( rho : double, p : double )
  return cmath.sqrt(gamma*p/rho)
end
