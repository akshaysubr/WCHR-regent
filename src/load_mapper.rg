-- Compile and link circuit.cc
local cwchr
do
  local root_dir = arg[0]:match(".*/") or "./"
  assert(os.getenv('LG_RT_DIR') ~= nil, "LG_RT_DIR must be given to build the mapper.")
  local runtime_dir = os.getenv('LG_RT_DIR') .. "/"
  local wchr_cc = root_dir .. "wchr.cc"
  local wchr_so
  if os.getenv('OBJNAME') then
    local out_dir = os.getenv('OBJNAME'):match('.*/') or './'
    wchr_so = out_dir .. "libwchr.so"
  elseif os.getenv('SAVEOBJ') == '1' then
    wchr_so = root_dir .. "libwchr.so"
  else
    wchr_so = os.tmpname() .. ".so" -- root_dir .. "wchr.so"
  end
  local cxx = os.getenv('CXX') or 'c++'

  local cxx_flags = os.getenv('CC_FLAGS') or ''
  cxx_flags = cxx_flags .. " -O2 -Wall -Werror"
  if os.execute('test "$(uname)" = Darwin') == 0 then
    cxx_flags =
      (cxx_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cxx_flags = cxx_flags .. " -shared -fPIC"
  end

  local cmd = (cxx .. " " .. cxx_flags .. " -I " .. runtime_dir .. " " ..
                 wchr_cc .. " -o " .. wchr_so)
  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. wchr_cc)
    assert(false)
  end
  terralib.linklibrary(wchr_so)
  cwchr = terralib.includec("wchr.h", {"-I", root_dir, "-I", runtime_dir})
end

return cwchr
