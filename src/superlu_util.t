import "regent"

local superlu = {}
do
  local superlu_library = "-lsuperlu"
  local superlu_include_dir = "/opt/SuperLU_5.2.1"
  local root_dir = arg[0]:match(".*/") or "./"
  local superlu_util_cc = root_dir .. "superlu_util.c"
  superlu_util_so = os.tmpname() .. ".so"
  local cc = os.getenv('CC') or 'cc'
  local cc_flags = "-O3 -Wall -Werror -std=c99"
  cc_flags = cc_flags .. " -I" .. superlu_include_dir
  local is_darwin = os.execute('test "$(uname)" = Darwin') == 0
  if is_darwin then
    cc_flags =
      (cc_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cc_flags = cc_flags .. " -shared -fPIC"
  end
  cc_flags = cc_flags .. " -lm -lblas " .. superlu_library 

  local cmd = (cc .. " " .. cc_flags .. " " .. superlu_util_cc .. " -o " .. superlu_util_so)
  print(cmd)

  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. superlu_util_cc)
    assert(false)
  end
  terralib.linklibrary(superlu_util_so)
  if is_darwin then
    terralib.linklibrary("libsuperlu.dylib")
  else
    terralib.linklibrary("libsuperlu.so")
  end
  superlu.c = terralib.includec("superlu_util.h", {"-I", root_dir, "-I", superlu_include_dir })
end

local c = regentlib.c
