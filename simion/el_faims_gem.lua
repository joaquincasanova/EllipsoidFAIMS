
function make_gem(ai,bi,ci,ao,bo,co,re)
	
  local filename = "el_faims.gem"
  local old_file_handle = io.output(filename)
	del = 4
	io.write ([[
pa_define(]] .. ao+del.. [[,]] .. bo+del.. [[, ]] .. co+del.. [[, planar, xzmirror)
  ]])
	io.write ([[
electrode(1) {fill{notin{sphere(0,0,0,]] .. ao .. [[,]] .. bo .. [[,]] .. co .. [[)}
				notin{cylinder(0,0,]] .. co+del .. [[,]] .. re .. [[,]] .. re .. [[,]] .. co+2*del.. [[)}}}
  ]])
	io.write ([[
electrode(2) {fill{within{sphere(0,0,0,]] .. ai .. [[,]] .. bi .. [[,]] .. ci .. [[)}}}
  ]])

  io.flush()  -- ensure file is immediately written to disk
  return filename
end

function test(ai,bi,ci,ao,bo,co,re)
  
  -- create GEM file of proper dimensions
  local gem_filename = make_gem(ai,bi,ci,ao,bo,co,re)

  -- convert GEM file to PA# file.
  local pasharp_filename = string.gsub(gem_filename, ".gem", ".pa#")
  simion.command("gem2pa " .. gem_filename .. " " .. pasharp_filename)

  -- refine PA# file.
  simion.command("refine " .. pasharp_filename)
end

-- Simulates all parameterizations.
function test_all()
  os.remove("*.pa*")
  os.remove("*.gem")
  test(127,127,127,147,147,147,5)
end

test_all()  -- do main function

--
