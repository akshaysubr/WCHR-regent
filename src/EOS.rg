import "regent"
local problem = require("problem")

local c     = regentlib.c
local cmath = terralib.includec("math.h")

-- Get gamma from the problem file
local gamma    = problem.gamma
local gamma_m1 = gamma - 1.0
local onebygm1 = 1.0 / gamma_m1

local Rgas     = problem.Rgas
local Cv       = Rgas * onebygm1
local Cp       = gamma * Cv

__demand(__inline)
task get_pressure( rho : double, rhoe : double )
  return gamma_m1 * rhoe
end

__demand(__inline)
task get_internal_energy( rho : double, p : double)
  return onebygm1 * p
end

__demand(__inline)
task get_temperature( rho : double, p : double )
  return p / (rho * Rgas)
end

__demand(__inline)
task get_sos( rho : double, p : double )
  return cmath.sqrt(gamma*p/rho)
end
