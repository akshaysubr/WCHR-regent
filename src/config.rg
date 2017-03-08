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
}

local cstring = terralib.includec("string.h")

terra print_usage_and_abort()
  c.printf("Usage: regent.py edge.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h               : Print the usage and exit.\n")
  c.printf("  -prefix {prefix} : Use {prefix} as prefix for file I/O.\n")
  c.printf("  -p {value}       : Set the number of parallel tasks to {value}.\n")
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
  self.parallelism = 4

  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    if cstring.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()
    elseif cstring.strcmp(args.argv[i], "-prefix") == 0 then
      i = i + 1
      cstring.strcpy(self.filename_prefix, args.argv[i])
      self.fileIO = true
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.parallelism = c.atoi(args.argv[i])
    end
    i = i + 1
  end

  self.prow, self.pcol = factorize(self.parallelism)
end

return Config
