import "regent"

local c       = regentlib.c
local cstring = terralib.includec("string.h")

require("fields")
require("SOE")
require("RHS")
require("boundary")

local problem = require("problem")
local Config  = require("config")

