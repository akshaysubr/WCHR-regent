import "regent"

local c = regentlib.c
local cmath = terralib.includec("math.h")

struct Config
{
  fileIO          : bool,
  filename_prefix : int8[256],
  parallelism     : int,
  prow            : int,
  pcol            : int,
  nstats          : int,
  restart         : bool,
  restart_count   : int,
}

local cstring = terralib.includec("string.h")

terra print_usage_and_abort()
  c.printf("Usage: regent.py main.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h                 : Print the usage and exit.\n")
  c.printf("  -prefix {prefix}   : Use {prefix} as prefix for file I/O.\n")
  c.printf("  -restart {vizdump} : Use {prefix} and restart from viz dump. (Should have same prow and pcol as the original run)\n")
  c.printf("  -stats {nstats}    : Print stats only every {nstats} steps.\n")
  c.printf("  -p {value}         : Set the number of parallel tasks to {value} (Default = 1).\n")
  c.printf("  -prow {value}      : [Optional] Set the number of parallel tasks in x decomposition to {value}.\n")
  c.printf("  -pcol {value}      : [Optional] Set the number of parallel tasks in z decomposition to {value}.\n")
  c.abort()
end

terra file_exists(filename : rawstring)
  var file = c.fopen(filename, "rb")
  if file == nil then return false end
  c.fclose(file)
  return true
end

terra factorize(parallelism : int)
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return { size_x, size_y }
end

terra Config:initialize_from_command( nx : int, ny : int, nz : int )
  self.fileIO = false
  cstring.strcpy(self.filename_prefix, "")
  self.parallelism = 1
  self.nstats = 1
  self.restart = false
  self.restart_count = 0

  var use_prow : bool = false
  var use_pcol : bool = false

  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    if cstring.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()
    elseif cstring.strcmp(args.argv[i], "-prefix") == 0 then
      i = i + 1
      cstring.strcpy(self.filename_prefix, args.argv[i])
      self.fileIO = true
    elseif cstring.strcmp(args.argv[i], "-restart") == 0 then
      i = i + 1
      self.restart = true
      self.restart_count = c.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-stats") == 0 then
      i = i + 1
      self.nstats = c.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.parallelism = c.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-prow") == 0 then
      i = i + 1
      self.prow = c.atoi(args.argv[i])
      use_prow = true
    elseif cstring.strcmp(args.argv[i], "-pcol") == 0 then
      i = i + 1
      self.pcol = c.atoi(args.argv[i])
      use_pcol = true
    end
    i = i + 1
  end

  regentlib.assert(use_prow == use_pcol, "Both prow and pcol should be specified or use -p {value} for autofactorization")
  if use_prow then
    self.parallelism = self.prow*self.pcol
  else
    self.prow, self.pcol = factorize(self.parallelism)
  end
end

return Config
