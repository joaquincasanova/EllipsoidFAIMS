-- SIMION Lua workbench user program for
-- Statistical Diffusion Simulation (SDS) Model.
--
-- REV-5a-2009-07-07
-- Runs under SIMION version 8.0.4 or above (as is).
--
-- This is a refined model of using gas kinetic numbers to simulate
-- ion mobility and diffusion via Stokes' Law, with diffusion derived
-- from kinetic jump statistics data, efficiently simulating pressures
-- in the atmospheric range.  This version supports multiple
-- ion definitions as well as temperature, pressure, and bulk gas
-- velocity fields (incorporated as array files or external functions).
--
-- See the accompanying README.html, collision_sds_documentation.pdf
-- files for full details on this program (as well as the papers
-- cited).
--
-- This program requires the existence of the following file:
--
--    mbmr.dat --  collision statistics file
--                 (see load_diffusion_statistics below)
--
-- Each file below is optional and ignored if omitted:
--
--    m_defs.dat -- ion parameter definition file
--                  holds ion diameter(nm), Ko(10-4m2V-1s-1) by ion mass
--                  (see load_mass_data below)
--
--    field definition files...Note: each file must have the same
--    dimensions as the associated potential array (see load_array below):
--
--    p_defs.dat    holds pressure field data (torr)
--    t_defs.dat    holds Temperature field data (K)
--    vx_defs.dat   holds Vx bulk gas velocity (m/s)
--    vy_defs.dat   holds Vy bulk gas velocity (m/s)
--    vz_defs.dat   holds Vz bulk gas velocity (m/s)
--
-- WARNING: If magnetic and electrostatic PAs overlap, then this
-- program must be active only in the magnetic one.  You can disable
-- SDS on a certain instance (e.g. #2) by doing this:
--
--   local SDS = simion.import 'collision_sds.lua'
--   SDS.instances[2] = false
--
-- Trajectory quality of 0 or negative is recommended for speed.
--
-- SOURCE:
--
--   This code and corresponding documentation are based on the original
--   SDS model by Dave Dahl, 2004-09-27, included in SIMION PRG format in
--   the supplementary material in the paper:
--
--   Anthony D. Appelhans and David A. Dahl.
--     SIMION ion optics simulations at atmospheric pressure.
--     International Journal of Mass Spectrometry.
--     Volume 244, Issue 1, 15 June 2005, Pages 1-14.
--     http://dx.doi.org/10.1016/j.ijms.2005.03.010
--
--   The SIMION PRG was converted to SIMION Lua by D.Manura, 2007-05,
--   and incorporated into SIMION 8 by permission of Appelhans and Dahl.
--   See the README.html for changes made.
--
--   Version 5 Supports capabilities for FAIMS simulations.
--
--   The current form is
--   (c) 2007-2009 Scientific Instrument Services, Inc. (Licensed
--   under SIMION 8.0)
--
simion.workbench_program()

-- detect program errors ( http://simion.com/issue/490 )
if checkglobals then checkglobals() end

-- check SIMION version.
assert(type(simion.VERSION) == 'table' and
       simion.VERSION >= simion.VERSION '8.0.4' or
       tostring(simion.VERSION) >= '8.0.4',
       "version 8.0.4 or above required")

local opt = ...

local M = {}
M.segment = {}


--##
--## SECTION: Common Constants
--##


-- Copy of function references (can increase performance and conciseness).
local abs   = math.abs
local min   = math.min
local max   = math.max
local sqrt  = math.sqrt
local log10 = math.log10
local exp   = math.exp
local modf  = math.modf
local floor = math.floor
local rand  = rand
local print = print
local assert= assert
local type  = type
local error = error
local pairs = pairs
local ipairs= ipairs


-- Physical constants used.
local ELEMENTARY_CHARGE = 1.602176462e-19  -- elementary charge (C)
local K_BOLTZMANN  = 1.3806503e-23         -- Boltzmann constant (JK-1)
local N_AVOGADRO   = 6.02214199e23         -- Avogadro's number
local MOL_VOLUME   = 22.413996e-3          -- Volume (m3) of one mol
                                           -- of ideal gas at 0 C, 101.325 kPa
local AMU_TO_KG    = 1.66053873e-27        -- mass of one u in kg
local PI           = math.pi               -- value of PI
local MM_HG_TO_PA  = 133.322               -- conversion of mmHg to Pa
local STP_TEMP     = 273.15                -- standard temperature (K)
local M_AIR        = 28.94515              -- Effective mass of air (u)
local D_AIR        = 0.366                 -- Effective diameter of air (nm)

--##
--## SECTION: Util
--##



--##
--## SECTION: I/O Utilities - general functions for input/output
--##


-- Returns (true/false) whether the file of the given name exists
-- (i.e. is readable at least).
local function file_exists(filename)
  local fh = io.open(filename)
  local is_exist = (fh ~= nil)
  if fh then fh:close() end
  return is_exist
end


-- Reads numeric data file with name filename into array that is
-- returned.  Commas and white-space are ignored.  Comments (which
-- begin with a semicolon and continue to the end of the line) are
-- ignored.  The table has an additional field 'ncols' indicating
-- the number of columns on the last non-empty line.
-- Raises on error.
local function read_file_numbers(filename)
  local ncols_cur = 0
  local ncols
  local function end_of_line()
    ncols = ncols_cur > 0 and ncols_cur or ncols
    ncols_cur = 0
  end

  local array = {}
  local f = assert(io.open(filename))
  while true do
    local v = f:read(1)
    if v then
      if v == ";" then -- comment line
        f:read() -- skip line
        end_of_line()
      elseif v == "," or v == " " or v == "\t" then -- delimiter
        -- ignore
      elseif v == "\n" or v == "\r" then
        end_of_line()
      else
        f:seek("cur", -1)
        local v2 = f:read("*number")
        if v2 then
          table.insert(array, v2)
          ncols_cur = ncols_cur + 1
        else
          f:close()
          error(v)
        end
      end
    else end_of_line() break end -- end of input
  end
  array.ncols = ncols
  f:close()
  return array
end


-- Same as read_file_numbers, except it just returns nil
-- (not raises an error) if file doesn't exist.
local function opt_read_file_numbers(filename)
  return file_exists(filename) and read_file_numbers(filename) or nil
end


-- Formatted print.
local function printf(...) print(string.format(...)) end


--##
--## SECTION: array - functions concerning arrays representing T,P,v.
--##



-- Creates function for interpolating given array value at point
-- (px,py,pz).
-- Function is of the form
--   array, px,py,pz -> value
--
-- This is a helper function for load_array
local function make_array_interp(array)
  local xsize, ysize, zsize = array.xsize, array.ysize, array.zsize
  local nx, ny = array.nx, array.ny
  local is3d = array.symmetry == '3dplanar[xyz]'

  -- Compute array offsets.
  -- Array data starts on 4th element
  local add = 4
  local offset1 = add                -- [0,0,0]
  local offset2 = add + 1            -- [1,0,0]
  local offset3 = add +     array.nx -- [0,1,0]
  local offset4 = add + 1 + array.nx -- [1,1,0]
  local add = array.nx * array.ny
  local offset5 = offset1 + add      -- [0,0,1]
  local offset6 = offset2 + add      -- [1,0,1]
  local offset7 = offset3 + add      -- [0,1,1]
  local offset8 = offset4 + add      -- [1,1,1]

  return function(array, px,py,pz)
    -- keep inside bounds
    if px >= xsize then
      if px > ysize then error" ERROR: ion's x coord outside array" end
      px = px - 0.000001
    end
    if py >= ysize then
      if py > ysize then error" ERROR: ion's y coord outside array" end
      py = py - 0.000001
    end
    if is3d and pz >= zsize then
      if pz > zsize then error" ERROR: ion's z coord outside array" end
      pz = pz - 0.000001
    end

    -- Compute closest lower-left-corner (llc) array indices and
    -- fraction parts in 2D
    local ipx, fx1 = modf(px)
    local ipy, fy1 = modf(py)
    local fx0, fy0 = 1-fx1, 1-fy1
    local index = ipy * nx + ipx

    -- Compute point weightings.
    local w1 = fx0 * fy0   -- [0,0]
    local w2 = fx1 * fy0   -- [1,0]
    local w3 = fx0 * fy1   -- [0,1]
    local w4 = fx1 * fy1   -- [1,1]
    if is3d then
      -- Compute closest lower-left-corner (llc) array indices and
      -- fraction parts in 3D
      local ipz, fz1 = modf(pz)
      local fz0 = 1-fz1
      index = index + ipz * (ny*nx)

      local w5 = w1 * fz1  -- [0,0,1]
      local w6 = w2 * fz1  -- [1,0,1]
      local w7 = w3 * fz1  -- [0,1,1]
      local w8 = w4 * fz1  -- [1,1,1]
            w1 = w1 * fz0  -- [0,0,0]
            w2 = w2 * fz0  -- [1,0,0]
            w3 = w3 * fz0  -- [0,1,0]
            w4 = w4 * fz0  -- [1,1,0]
      return
        array[index + offset1] * w1 +
        array[index + offset2] * w2 +
        array[index + offset3] * w3 +
        array[index + offset4] * w4 +
        array[index + offset5] * w5 +
        array[index + offset6] * w6 +
        array[index + offset7] * w7 +
        array[index + offset8] * w8
    else -- 2D
      return
        array[index + offset1] * w1 +
        array[index + offset2] * w2 +
        array[index + offset3] * w3 +
        array[index + offset4] * w4
    end
  end
end


-- Loads 2D or 3D array of data from text file with name filename.
-- For example, this can be used as a pressure, temperature, or
-- velocity map (to allow these parameters to vary as a function of
-- position).
--
-- The file must have the format as described in read_file_numbers and
-- contain these values:
--
--   nx,ny,nz -- first three file entries are array dimensions
--            --   (number of points in x, y, and z directions)
--   val      -- rest of array is linear list of values
--   val      --   (suggest one per line), starting at
--   val      --   x,y,z=0,0,0 and scanning in x, tehen by y, then by z
--   ...
--
-- Note that
--
--   nx = 0    -- signifies an empty array (array is ignored)
--   nz = 0    -- 2D cylindrical array
--   nz = 1    -- 2D planar array
--
-- The result is a table containing these fields:
--
--   symmetry -- '2dplanar[xy]', '2dcylindrical[x]', or '3dplanar[xyz]'
--   nx       -- number of points in x dimension
--   ny       -- number of points in y dimension
--   nz       -- number of points in z dimension
--   xsize    -- x size (nx - 1)
--   ysize    -- y size (ny - 1)
--   zsize    -- z size (nz - 1)
--
-- Returns nil if the file does not exist or array is empty.
--
local function load_array(filename)
  local array = opt_read_file_numbers(filename)
  if array and array[1] == 0 then array = nil end

  if array then
    -- Get dimensions.
    local x_dim,y_dim,z_dim = array[1],array[2],array[3]
    printf("SDS Note: Loading array %s: nx=%d, ny=%d, nz=%d",
           filename, x_dim, y_dim, z_dim)
    array.nx = array[1]
    array.ny = array[2]
    array.nz = array[3]
    array.symmetry = (array.nz > 1) and '3dplanar[xyz]' or '2dplanar[xy]'
    if array.nx < 1 then error("X dimension invalid") end
    if array.ny < 1 then error("Y dimension invalid") end
    if array.nz == 0 then   -- cylindrical
      array.symmetry = '2dcylindrical'
      array.nz = 1
      print("SDS Note: 2D cylindrical arrays")
    end
    if array.nz < 1 then error("Z dimension invalid") end

    -- Compute limits and borders.
    array.xsize = abs(array.nx - 1)  -- max addressable x index (0 - nx-1)
    array.ysize = abs(array.ny - 1)  -- max addressable y index (0 - ny-1)
    array.zsize = abs(array.nz - 1)  -- max addressable z index (0 - nz-1)

    -- Allow object to be callable, returning interpolated values.
    setmetatable(array, {__call = make_array_interp(array)})
  end

  return array
end
M.load_array = load_array


--##
--## SECTION: Diffusion Statistics
--##


-- Number of collisions represented in distribution data (see
-- load_diffusion_statistics).
local N_DIST_COLLISIONS = 100000

-- Number of ICDFs in represented in distribution data (see
-- load_diffusion_statistics)
local N_DIST = 5

-- Number of data points per ICDF.
-- (Note: the extra point, 1002, might be unnecesary.)
local N_DIST_POINTS = 1002


-- Loads and returns diffusion statistics data from file.
--
-- This data represents a number of inverse cumulative probability
-- density functions (ICDFs), one for each of a series of mass ratios
-- (mass_ion/mass_gas), which are powers of 10 from 1 to 10,000
-- (referred to as functions #1 - #5).
--
-- Each function represents the distance that a particle of the given
-- mass ratio travels (between starting and ending points) upon doing
-- N_DIST_COLLISIONS collisions, assuming a mean-free-path of 1 unit.
-- The function is represented by N_DIST_POINTS data points to cover
-- percentiles 0% to 100% in increments of 0.1%.
--
-- The data are fully randomized Maxwell-Boltzmann (MB) statistics:
--   MB randomized initial energy (normalized to average velocity of 1).
--   Randomized initial direction.
--   Poisson randomized collision distance (normalize to average of 1).
--   MB randomized collision gas energies (normalized to have same mean
--     energy of the ions)
--   Randomized equally probable impact points.
--
-- N_DIST_COLLISIONS is also stored in the returned table.
--
-- Raises on error.
local function load_diffusion_statistics(filename)
  local raw = read_file_numbers(filename)  -- raises
  assert(#raw == N_DIST_POINTS * N_DIST)

  -- Split ICDFs into separate arrays.
  local stats = {}
  for i=1,N_DIST do
    local icdf = {}
    for j=1,N_DIST_POINTS do
      icdf[j] = raw[(i-1) * N_DIST_POINTS + (j-1) + 1]
    end
    stats[i] = icdf
  end
  stats.N_DIST_COLLISIONS = N_DIST_COLLISIONS

  return stats
end


-- Compute a random walk jump distance for diffusion effect.  This is
-- the distance traveled (between starting and ending points) by a
-- particle after having done N_DIST_COLLISIONS collisions, assuming
-- log10(mass_ion/mass_gas) = log_mr_ratio and normalized conditions
-- described in load_diffusion_statistics.
--
-- This value is interpolated from statistics in stats, which
-- represents inverse cumulative density function (ICDF).  It involves
-- a log-log scale interpolation of ICDFs at the nearest mass ratios
-- (the value of each ICDFs is itself interpolated as well).
local function diff_dist_steps(stats, log_mr_ratio)
  -- Select uniformly random percentile * 10 (0, 1000) for the random
  -- jump distance.
  local n = rand() * (N_DIST_POINTS - 2)

  -- Select the two ICDFs to interpolate between.  We'll interpolate
  -- between ICDFs icdf1 and icdf2.  These are selected based on mass
  -- ratio.
  local iicdf = (log_mr_ratio <= 1) and 1 or
                (log_mr_ratio <= 2) and 2 or
                (log_mr_ratio <= 3) and 3 or 4
  local icdf1,icdf2 = stats[iicdf], stats[iicdf+1]

  -- Select the indices in each ICDF to interpolate between.
  local ilow = floor(n) + 1; local ihigh = ilow + 1

  -- Interpolate inside both ICDFs.
  local weight = n % 1
  local dist1_steps = (icdf1[ihigh] - icdf1[ilow]) * weight + icdf1[ilow]
  local dist2_steps = (icdf2[ihigh] - icdf2[ilow]) * weight + icdf2[ilow]

  -- Interpolate between both ICDFs (log-log interpolation).
  dist1_steps, dist2_steps = log10(dist1_steps), log10(dist2_steps)
  local weight = log_mr_ratio - (iicdf - 1)
  local dist_steps = 10^((dist2_steps - dist1_steps) * weight + dist1_steps)

  return dist_steps
end


--##
--## SECTION: Mass Data (massdata)
--##


-- Gets scaling constant air_to_gas such that Kogas = Koair *
-- air_to_gas.  Given diameter d_ion (nm) and mass mass_ion (u)
-- of ion and given diameter d_gas (nm) and mass mass_gas (u) of gas.
local function get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)
  -- Reduced mass with respect to air and gas
  local reduced_mass_air = mass_ion * M_AIR    / (mass_ion + M_AIR)
  local reduced_mass_gas = mass_ion * mass_gas / (mass_ion + mass_gas)

  -- Scaling constant to convert from Koair to Kogas.
  local air_to_gas =
    ((d_ion + D_AIR) / (d_ion + d_gas) )^2 *
    sqrt(reduced_mass_air / reduced_mass_gas)

  return air_to_gas
end


-- Estimate ion diameter d_ion (nm) from mass mass_ion (u) and
-- (optionally) mobility (ko).  Given diameter d_gas (nm) and mass
-- mass_gas (u) of gas.
--
local function estimate_d_ion(mass_ion,ko, d_gas,mass_gas)
  -- Get rough estimate of ion diameter (in nm) from mass_ion (u).
  -- Emperical formula noted in paper.
  local d_ion = 0.120415405 * mass_ion^(1/3)

  if ko then -- refinement (if ko given)
    -- Compute Koair by scaling Ko(gas).
    local Koair = ko / get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)

    -- Estimate d (nm) from Ko (10-4 m2 V-1 s-1) from ko.
    local C = 1.0e5  -- converts to (10-9 m2 V-1 s-1)
    local logkm = log10(Koair * C)
    -- Emperical formula noted in the paper.
    d_ion = 10^(3.0367 - 0.8504 * logkm + 0.1137 * logkm^2 - 0.0135 * logkm^3)
  end

  return d_ion
end


-- Estimate ion mobility (ko) from mass mass_ion (u) and diameter
-- d_ion (nm) of ion.  Given diameter d_gas (nm) and mass mass_gas
-- (u) of gas.
--
local function estimate_ko(mass_ion,d_ion, d_gas,mass_gas)
  -- Estimate Koair (10-4 m2 V-1 s-1) from d_ion
  -- Emperical formula noted in paper.
  local logdm = log10(d_ion)
  local Koair =
    1.0e-5 * 10^(4.9137 - 1.4491*logdm - 0.2772*logdm^2 + 0.0717*logdm^3)

  -- Compute Kogas by scaling Koair (10-4 m2 V-1 s-1)
  local ko = Koair * get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)

  return ko
end


-- Adds ions STP average velocity Vo and mean free path MFPo data to
-- given mass record mdata.  The data is computed from the other data
-- in the mass record.  Given diameter d_gas (nm) and mass mass_gas
-- (u) of gas.
-- returns mfpmass, vmass (the Vo and MFPo values added).
local function massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  local mass_ion, d_ion, ko = mdata.mass, mdata.d, mdata.ko

  -- Mean ion speed at STP (mm/usec)
  local MM_USEC__M_S = 1.0e-3  -- (mm/usec)/(m/s)
  local Vk = sqrt(8 * K_BOLTZMANN * STP_TEMP / PI / AMU_TO_KG) * MM_USEC__M_S
  -- print('DEBUG:check: ', Vk, '=', 2.404850043)

  -- Thermal velocity of ion (mm/usec)
  local Vio = Vk * sqrt(1 / mass_ion)

  -- Thermal velocity of gas molecule (mm/usec)
  local Vgo = Vk * sqrt(1 / mass_gas)

  -- Particles per volume, No, (n/mm3) at STP.
  local MM3_PER_M3 = 1.0e9     -- (mm3/m3)
  local No = N_AVOGADRO        -- (n/mol)
           / MOL_VOLUME        -- (m3/mol)
           / MM3_PER_M3

  -- Collision frequency (collisions/usec)
  -- Collision frequency constant Fk 
  local MM2_PER_NM2 = 1.0e-12  -- (mm2/nm2)
  local Fio = MM2_PER_NM2 *  No * PI *
    ((sqrt(2)-1/4) * ((d_gas + d_ion)/2)^2 * Vio + (1/4) * d_ion^2 * Vgo)

  -- Ion's mean free path (mm) at STP
  local Lio = Vio / Fio

  local vmass   = Vio      -- ion's STP mean velocity (mm/usec), Vo
  local mfpmass = Lio      -- ion's STP mean free path, MFPo

  -- Store Vo and MFPo in mass record.
  mdata.vo   = vmass
  mdata.mfpo = mfpmass

  return mfpmass, vmass
end


-- Gets data (MFPo, Vo, and Ko, alpha, beta) for given mass mass_ion
-- (u) from massdata.  Fills in mass record with estimated values if
-- no match.
local function massdata_for_mass(massdata, mass_ion)
  -- Return any record that matches.
  for i,mdata in ipairs(massdata) do
    if mass_ion == mdata.mass then
      return mdata.mfpo, mdata.vo, mdata.ko,
             mdata.alpha, mdata.beta  -- FAIMS
    end
  end
  -- Otherwise, estimate values...

  local d_gas, mass_gas = massdata.d_gas, massdata.mass_gas

  -- Estimate ion diameter (nm) from mass_ion
  local d_ion = estimate_d_ion(mass_ion,nil, d_gas,mass_gas)

  -- Estimate Ko (10-4 m2 V-1 s-1) from mass_ion and d_ion.
  local ko = estimate_ko(mass_ion,d_ion, d_gas,mass_gas)

  -- BEGIN FAIMS
  -- Estimate mobility factors alpha and beta
  local alpha, beta = 0, 0
  -- END FAIMS

  -- Create new mass record.
  local mdata = {
    mass       = mass_ion,
    d          = d_ion,
    ko         = ko,
    -- BEGIN FAIMS
    alpha      = alpha,
    beta       = beta,
    -- END FAIMS
    estimation = 'dk'
  }
  -- add Vo and MFPo also to mdata
  local mfpmass, vmass = massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  massdata[#massdata+1] = mdata

  return mfpmass, vmass, ko,
         alpha, beta  -- FAIMS
end



-- Fill in missing data in massdata mass records, after diameter d_gas
-- (nm) and mass mass_gas (u) of background gas is known.
local function massdata_complete_records(massdata, d_gas,mass_gas)
  -- Preserve.
  massdata.d_gas, massdata.mass_gas = d_gas,mass_gas

  -- Update each record.
  for i,mdata in ipairs(massdata) do
    local mass_ion, d_ion, ko = mdata.mass, mdata.d, mdata.ko

    -- Compute d or Ko if necessary.
    -- Update mass record.
    if ko * d_ion ~= 0 then     -- both defined
      mdata.estimation = ''
    elseif d_ion == 0 and ko == 0 then -- both undefined
      mdata.d = estimate_d_ion(mass_ion,nil, d_gas,mass_gas)  -- ko = nil
      mdata.ko = estimate_ko(mass_ion,mdata.d, d_gas,mass_gas)
      mdata.estimation = 'dk'
    elseif d_ion == 0 then      -- d_ion undefined
      mdata.d = estimate_d_ion(mass_ion,ko, d_gas,mass_gas)
      mdata.estimation = 'd'
    elseif ko == 0 then         -- Ko undefined
      mdata.ko = estimate_ko(mass_ion,d_ion, d_gas,mass_gas)
      mdata.estimation = 'k'
    end

    -- Add Vo and Mean free path to record.
    massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  end
end


-- Print data for all mass records.
local function massdata_print(massdata, tlocal,plocal)
  -- Print mass data.

  -- sort data by increasing mass
  local mass_idxs = {}
  for i=1,#massdata do mass_idxs[i] = i end
  table.sort(mass_idxs,
    function(a,b) return massdata[a].mass < massdata[b].mass end)

  -- Print mass data.
  print()
  print("SDS Mass definitions loaded or created (sorted by mass)")
  print("  N , mass, dia,  Ko Mobility  ,alpha, beta, MFPo, Avg Vo")
  print(" (n),(amu),(nm),(10-4 m2V-1s-1),(Td-2),(Td-4),(mm),(mm/usec)")
  local missing = false
  for _,i in ipairs(mass_idxs) do
    local mdata = massdata[i]
    -- add line to new mass data output
    printf(
      " n=%d, m=%g, d=%g(%s), Ko=%g(%s), a=%g, b=%g, MFPo=%g, Vo=%g %s", 
      i,
      mdata.mass,
      mdata.d,
      (mdata.estimation:find'd' and "est" or "def"),
      mdata.ko,
      (mdata.estimation:find'k' and "est" or "def"),
      mdata.alpha or 0,
      mdata.beta or 0,
      mdata.mfpo,
      mdata.vo,
      mdata.estimation ~= '' and '[*]' or ''
    )
    missing = missing or mdata.estimation ~= ''
  end
  if missing then
    print(" [*] WARNING: Estimations made due to absent data in m_defs.dat")
  end

  -- Print mass data at local conditions.
  printf("")
  printf("SDS Local values for ions (sorted by mass) at user defined")
  printf("   Pressure = %g Torr, Temperature = %g K", plocal, tlocal)
  printf(" mass ,  K Mobility   , MFP,  Avg V  ,B1 MFP/rB1=1")
  printf(" (amu),(10-4 m2V-1s-1),(mm),(mm/usec),(gauss)")
  for _,i in ipairs(mass_idxs) do
    local mdata = massdata[i]

    -- Correct for local temperature and pressure.
    local vlocal   = mdata.vo * sqrt(tlocal / STP_TEMP)
    local mfplocal = mdata.mfpo * (760 / plocal) * (tlocal / STP_TEMP)
    local klocal   = mdata.ko   * (760 / plocal) * (tlocal / STP_TEMP)

    -- Magnetic field B when MFP = r_cyclotron (thermal).
    local B1 = 1439.74 *
               sqrt(mdata.mass * speed_to_ke(vlocal, mdata.mass)) / mfplocal

    printf(" m=%g, K=%g, MFP=%g, v=%g, B1=%g",
           mdata.mass, klocal, mfplocal, vlocal, B1)
  end

  print()
end


-- Load mass definitions (ion diameters and mobility constants) keyed to mass
-- from file with name filename.
--
-- The file must have the format as described in read_file_numbers and
-- contain tuples of values as follows:
--
--   mass(u),d(nm),Ko(10-4 m2 V-1 s-1),alpha(Td-2),beta(Td-4)
--   mass(u),d(nm),Ko(10-4 m2 V-1 s-1),alpha(Td-2),beta(Td-4)
--   ...
--
-- This routine creates an array of mass records.  Records are tables
-- with these fields:
--
--   mass  -- ion mass (u)
--   d     -- ion diameter d (nm)
--   ko    -- ion reduced mobility Ko (10-4 m2 V-1 s-1) (at STP).
--   alpha -- (Td-2)  (optional)
--   beta  -- (Td-4)  (optional)
--
-- The alpha and beta constants allow mobility to vary as a function
-- to field strength (e.g. FAIMS).  Set them to zero to prevent this.
-- Mobility at a given E/N (electric field per number density of
-- particles) is approximated as
--
--   K = ko * (1 + alpha * (E/N)^2 + beta * (E/N)^4 + ...)
--
-- Records can also contain these fields (after calling
-- massdata_complete_records):
--
--   vo       -- ion velocity Vo (mm/usec)
--   mfpo     -- ion MFPo (mm)
--   estimation -- which of ion d and Ko are estimated rather than defined:
--                 '', 'd', 'k', 'dk'
--
-- The containing array can also contain these fields (after calling
-- massdata_complete_records):
--
--   d_gas    -- background gas diameter (nm)
--   mass_gas -- background gas mass (u)
--
-- returns empty list if file not exist.
--
local function load_massdata(filename)
  local raw = opt_read_file_numbers(filename)
  local massdata = {}
  -- Check mass records.
  if raw then
    local ncols = raw.ncols --FIX
    if (ncols ~= 3 or ncols ~= 5) and #raw % ncols ~= 0 then
      error(filename .. " does not contain three or five numbers per line.")
    end
    local ir, i = 1, 1
    while raw[ir] and raw[ir] ~= 0 do
      massdata[i] = {
        mass  = abs(raw[ir]),
        d     = abs(raw[ir + 1]),
        ko    = abs(raw[ir + 2]),
        alpha = ncols >= 5 and raw[ir + 3] or nil,
        beta  = ncols >= 5 and raw[ir + 4] or nil
      }
      ir = ir + ncols
      i = i + 1
    end
    -- check for duplicated masses
    local found = {}
    for _,mdata in ipairs(massdata) do
      local mass = mdata.mass
      if found[mass] then error("SDS ERROR: Mass duplicated: #", mass) end
      found[mass] = true
    end
  else
    printf("SDS Warning: mass definitions \"%s\" not found", filename)
  end
  return massdata
end


--##
--## SECTION: More Locals
--##


-- Collision statistics data.
local s_dist = load_diffusion_statistics("mbmr.dat")

-- Data for diffusion and mobility parameters keyed to mass.
local massdata = load_massdata("m_defs.dat")

-- User supplied initialization function.
-- If not nil, this function is called at the beginning of SDS initialization.
-- A typical application of this is to intialize SDS adjustable variables
-- based on other adjustable variables. Example:
--
--   adjustable pressure_pascals = 101325
--   function SDS.init()
--     adjustable SDS_pressure_torr = pressure_pascals * (760/101325)
--   end
local init = nil

-- Definitions for pressure, temperature, and background gas velocity.
-- These can be functions (e.g. for an analytical equation) or a data
-- array (e.g. for interpolation from grid data) or nil (not used).
--
-- pressure function (torr)
--
local pres_defs = load_array("p_defs.dat")
--
-- temperature function (deg K)
--
local temp_defs = load_array("t_defs.dat")
--
-- vx bulk gas velocity array (m/s).
-- Note: vx,vy,vz, if defined via an array, are oriented with respect
-- to the PA (not the workbench).
--
local velx_defs = load_array("vx_defs.dat")
--
-- vy bulk gas velocity array (m/s)
--
local vely_defs = load_array("vy_defs.dat")
--
-- vz bulk gas velocity array (m/s)
--
local velz_defs = load_array("vz_defs.dat")
--
-- bulk gas velocity function (m/s) for all three components
-- (vx,vy,vz).  if set to not nil, this overrides the above three
-- values.
-- Note: vx,vy,vz, if defined via a function, are oriented with respect
-- to the workbench (not the PA), which is different from the above.
--
local vel_defs
--
-- The coordinate system orientation in which the velocity vector
-- returned by the above velocity objects is understood.
--
--  "mm"     -- workbench coordinates
--  "gu"     -- array volume coordinates
--  "abs_gu" -- absolute array coordinates (e.g removing symmetry)
--
-- (See section L.1.6 "Units and Coordinate Conventions" of
-- the printed SIMION manual for details.)
local vel_coords
--
-- The coordinate system in which pressure function is evaluated.
local pres_coords
-- The coordinate system in which temperature function is evaluated.
local temp_coords

-- will be set to true if any local parameters above are defined.
local is_local_params = false

-- will be set to one of the arrays to use a template for dimensions.
local array_template


---- Data for each ion.
---- For each variable t below, t[i] is the data for ion number i.
---- Note that STP standard for "standard temperature and pressure".
--
-- Ion's Stokes' law damping (usec-1) at STP.
-- This is converted from ion mobility.
--
local ions_STP_damping = {}
--
-- Ion's mean free path (mm) at STP.
--
local ions_STP_mfp_mm = {}
--
-- Ion's average thermal velocity (mm/usec) at STP.
--
local ions_STP_Vo_mm_per_usec = {}
--
-- Ion's Stoke's law damping (usec-1) at the local temperature and
-- pressure.  This is calculated in part from ions_STP_damping.
--
local ions_local_damping = {}
--
-- Ion's local pressure, temp corrected mean free path (mm).
-- This is calculated in part from ions_STP_mfp_mm.
--
local ions_local_mfp_mm = {}
--
-- Ion's average thermal velocity rate (mm/usec) at local temperature.
-- This is calculated in part from ions_STP_Vo_mm_per_usec.
--
local ions_local_V_mm_per_usec = {}
--
-- Ion's log10(mass_ion/mass_gas), where mass_ion/mass_gas is ratio of
-- ion mass to gas mass.
--
local ions_log_mr_ratio = {}

-- BEGIN FAIMS
--
-- Ion's alpha mobility constant (Td-2).
-- K(E/N) / K_0 = 1 + alpha (E/N)^2 + beta (E/N)^4 + ...
--
local ions_alpha = {}
--
-- Ion's beta mobility constant (Td-4).
-- K(E/N) / K_0 = 1 + alpha (E/N)^2 + beta (E/N)^4 + ...
--
local ions_beta = {}
--
-- ions_alphabeta[n] true iff ions_alpha[n] ~= 0 or ions_beta[n] ~= 0
--
local ions_alphabeta = {}
-- END FAIMS

-- System parameters at current ion position.
-- These apply to the current ion.
local local_pressure_torr   -- local pressure (torr)
local local_temperature_K   -- local temperature (K)
local local_vx_mm_per_usec  -- local bulk gas velocity components (mm/usec)
local local_vy_mm_per_usec  -- ..
local local_vz_mm_per_usec  -- ..
-- BEGIN FAIMS
local local_dvoltsx_Vmm     -- local field (V/mm)
local local_dvoltsy_Vmm
local local_dvoltsz_Vmm
local local_dvoltsx_last_Vmm = {}   -- local field (V/mm) at end of time step
local local_dvoltsy_last_Vmm = {}   --   keyed by ion number.
local local_dvoltsz_last_Vmm = {}
-- END FAIMS


-- The user program is enabled when instances[ion_instance] = true.
-- Set to false to disable SDS on certain PA instance numbers
-- (e.g. when electrostatic and magnetic PAs overlap, you only want
-- SDS to be active only on the magnetic PA).
local instances = {}
for i=1,200 do instances[i] = true end -- all instances enabled by default.


if sim_trajectory_quality > 0 then
  print("SDS Warning: Trajectory quality (TQual) <= 0 is " ..
        "recommended for faster speed.")
end


--##
--## SECTION: SDS Util
--##


-- Checks that array has same dimension as previous array.  Also sets
-- is_local_params and array_template.  This is a helper function for
-- check_all_arrays.
local function check_array(array)
  if type(array) == 'userdata' then -- possibly SIMION PA
    assert(type(simion.VERSION) == 'table' and
           simion.VERSION >= simion.VERSION '8.0.5-TEST19',
           "version 8.0.5-TEST19 or above required")
  end

  if array then is_local_params = true end
  if not array_template and type(array) == 'table' then
    array_template = array
  end
  if type(array) ~= 'table' then return end
  if array.nx ~= array_template.nx then
     error("ERROR: nx dimension of arrays don't match") end
  if array.ny ~= array_template.ny then
     error("ERROR: ny dimension of arrays don't match") end
  if array.nz ~= array_template.nz then
     error("ERROR: nz dimension of arrays don't match") end
end


-- checks data arrays: check and init for P,T, Vel, and massdata data arrays.
-- Has side-effects.
local function check_all_arrays()
  -- Check all arrays (these can be nil)
  local arrays = {pres_defs, temp_defs, velx_defs,
                  vely_defs, velz_defs, vel_defs}
  for _,array in pairs(arrays) do   -- note: careful with nil's in table.
    check_array(array)
  end

  assert(type(velx_defs) ~= 'function')
  assert(type(vely_defs) ~= 'function')
  assert(type(velz_defs) ~= 'function')
end


-- Update local data of ion number nion for local temperature and
-- pressure (i.e. ions_local_* variables).  This uses system local
-- variables (local_temperature_K, local_pressure_torr,
-- local_dvolts_Vmm) as well as ion's STP variables (ions_STP_mfp_mm
-- and ions_STP_Vo_mm_per_usec).
local function update_ions_local(nion)
  local mfpmass = ions_STP_mfp_mm[nion]
  local vmass   = ions_STP_Vo_mm_per_usec[nion]

  -- Correct for local temperature and pressure.
  local t_ratio = local_temperature_K / STP_TEMP          -- local T
  local pt_ratio = t_ratio * (760 / local_pressure_torr)  -- local P
  local o_ratio = 1  -- other corrections

  -- BEGIN FAIMS
  --
  -- E/N in units of Td = 1E-21 (V m-1)(m3)
  --
  -- Particles per volume, No, (n/mm3) at STP.
  -- WARNING: using ideal gas law (might not be the most accurate).
  local alpha, beta = ions_alpha[nion], ions_beta[nion]
  if alpha ~= 0 or beta ~= 0 then
    local MM3_PER_M3 = 1.0e9     -- (mm3/m3)
    local No = N_AVOGADRO        -- (n/mol)
             / MOL_VOLUME        -- (m3/mol)
    -- Particles per volumes, No, (n/mm3) at local conditions.
    local N = No / pt_ratio      -- (n/m3)
    -- E
    local fx = local_dvoltsx_Vmm
    local fy = local_dvoltsy_Vmm
    local fz = local_dvoltsz_Vmm
    local field_Vmm = sqrt(fx*fx+fy*fy+fz*fz) -- (V/mm)
    -- E/N ratio.
    local en = (1E3 * 1E21) * field_Vmm / N  -- Td = 1E-21 V m^2
    -- Ratio K(E/N) / K_0 = 1 + alpha*(E/N)^2 + beta*(E/N)^4 + ...
    local en2 = en*en; local en4 = en2*en2
    o_ratio = 1 + alpha*en2 + beta*en4
  end
  -- END FAIMS

  ions_local_damping[nion] = ions_STP_damping[nion] / pt_ratio / o_ratio
  ions_local_mfp_mm[nion] = mfpmass * pt_ratio
  ions_local_V_mm_per_usec[nion] = vmass * sqrt(t_ratio)
end



--##
--## SECTION: SIMION Adjustable Variables (used in segments).
--##


-- Whether SDS effects are enabled
-- 0=no
-- 1=yes
adjustable SDS_enable = 1

-- Whether diffusion effects enabled.
-- 0=no, 1=yes (default yes).
-- It can sometimes be helpful to disable this to observe the behavior
-- without diffusion effects, such as to observe the average expected
-- behavior using a single ion.
adjustable SDS_diffusion = 1

-- BEGIN FAIMS
--
-- Whether FAIMS mode is in effect (1=yes,0=no)
-- This enables optimizations for FAIMS mode.
-- Use this only for strong RF fields, where the mobility constant is
-- a function of the field intensity and the ions reach terminal velocity
-- within a very small fraction of the RF period time.
-- See our FAIMS paper for details and sections of the code that test
-- SDS_faims_mode.
-- 
-- References:
--   http://www.faims.com/
--   http://en.wikipedia.org/wiki/High-field_asymmetric_waveform_ion_mobility_spectrometry
adjustable SDS_faims_mode = 0
-- END FAIMS

-- Mass of background gas particle (u).  Default assumes normal air
-- mixture.
adjustable SDS_collision_gas_mass_amu = 28.94515

-- Effective diameter of background gas molecules (nm).  Default
-- assumes air.
adjustable SDS_collision_gas_diameter_nm = 0.366

-- Background gas pressure (torr).  Default is 760 Torr = 101,325 Pa
-- (atmospheric pressure).
-- This will be overriden by pres_defs (if defined)
adjustable SDS_pressure_torr = 760.00

-- Background gas temperature (K).  Default is standard temp (25 C).
-- This will be overriden by temp_defs (if defined)
adjustable SDS_temperature_K = 298.15

-- Background gas mean velocity in x (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vx_m_per_sec  = 0.0

-- Background gas mean velocity in y (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vy_m_per_sec  = 0.0

-- Background gas mean velocity in z (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vz_m_per_sec  = 0.0

-- Scaling factor to apply on pressure array fields (pres_defs) if
-- defined.  Default is 1 (no scaling).
adjustable SDS_P_field_scale_factor = 1.0

-- Scaling factor to apply on temperature array fields (temp_defs) if
-- defined.  Default is 1 (no scaling).
adjustable SDS_T_field_scale_factor = 1.0

-- Scaling factor to apply on velocity (Vx, Vy, Vz) array fields
-- (vel*_defs) if defined.  Default is 1 (no scaling).
adjustable SDS_v_field_scale_factor = 1.0

-- Force minimum time step in usecs.
-- 0 (the default) disables this. >0 forces code speedup.
adjustable SDS_min_time_step_usec = 0.0

-- The following options pertain to the optional "quick RF cycle" feature,
-- which can improve speed, especially if a large number of
-- RF cycles must be simulated over a short distance traveled.
-- It simulates a single RF cycle
-- normally (over a time SDS_quick_period) but without
-- diffusion and computes the ion displacement
-- over that cycle.  Then it multiplies that displacement by
-- (SDS_quick_cycles + 1) and applies diffusion (if enabled)
-- expected for a time period of (SDS_quick_cycles + 1) RF cycles.
-- See also README.html file.
--
-- Number of RF cycles to simulate quickly after each
-- RF cycle simulated normally.
-- Default is 0 (no quick RF cycles).
-- Set to -1 to automatically choose value (chosen in such a way to
-- advance particle no more than one grid unit from mobility effect).
adjustable SDS_quick_cycles = 0

-- Time period of quick RF cycle (in usec).
-- This is only used if SDS_quick_cycles is non-zero.
-- Normally this should be exactly the same as the RF period.
-- Its recommended that the RF waveform code set this automatically
-- to ensure consistency.
adjustable SDS_quick_period = 0

--##
--## SECTION: SIMION Segments (or Code Designed to Run in them)
--##


-- Print coordinates of ion and local pressure, temperature, and gas
-- velocity.  Useful for debugging.
--
-- This is designed to be called in certain SIMION segments, such as
-- initialize, other_actions, terminate, or accel_adjust.
local function print_coordinates()
  printf("-----")
  printf("  Ion at array coords (gu): x=%f, y=%f, z=%f",
         ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu)
  printf("  Ion at PA coords (gu): x=%f, y=%f, z=%f",
         ion_px_gu, ion_py_gu, ion_pz_gu)
  printf("  Ion at workbench coords (mm): x=%f, y=%f, z=%f",
         ion_px_mm, ion_py_mm, ion_pz_mm)
  printf("  Bulk gas vel (m/s): vx=%f, vy=%f, vz=%f",
         local_vx_mm_per_usec * 1.0e3, local_vy_mm_per_usec * 1.0e3,
         local_vz_mm_per_usec * 1.0e3)  -- scale to m/s
  printf("  Local P(torr)=%f, T(K)=%f",
         local_pressure_torr, local_temperature_K)
end


-- Print SDS parameters to log.
-- (useful to ensure desired parameters are active)
local function print_sds_parameters()
  print("SDS collision gas diameter (nm)=", SDS_collision_gas_diameter_nm)
  print("SDS collision gas mass (u)=", SDS_collision_gas_mass_amu)

  local s= "SDS pressure (torr)= "
  if pres_defs then
    s = s .. type(pres_defs)
    if SDS_P_field_scale_factor ~= 1 then
      s = s .. ", scaled by factor=" .. SDS_P_field_scale_factor
    end
  else
    s = s .. SDS_pressure_torr
  end
  print(s)

  local s = "SDS temperature (K)= "
  if temp_defs then
    s = s .. type(temp_defs)
    if SDS_T_field_scale_factor ~= 1 then
      s = s .. ", scaled by factor=" .. SDS_T_field_scale_factor
    end
  else
    s = s .. SDS_temperature_K
  end
  print(s)

  local s = "SDS velocity (m/s)= "
  if vel_defs or velx_defs or vely_defs or velz_defs then
    s = s ..
      (vel_defs  and " vxyz=" .. type(vel_defs)  or '') ..
      (velx_defs and " vx="   .. type(velx_defs) or '') ..
      (vely_defs and " vy="   .. type(vely_defs) or '') ..
      (velz_defs and " vz="   .. type(velz_defs) or '')
    if SDS_v_field_scale_factor ~= 1 then
      s = s .. ", scaled by factor=" .. SDS_v_field_scale_factor
    end
  else
    s = s .. SDS_vx_m_per_sec .. ', ' ..  SDS_vy_m_per_sec .. ', '
          .. SDS_vz_m_per_sec
  end
  print(s)

  if not is_local_params then
    print("SDS Note: No local P, T, or v field defined.")
  end

  if SDS_min_time_step_usec ~= 0 then
    print("SDS min time step (usec)=", SDS_min_time_step_usec)
  end
end


-- Rotate velocity vector (x,y,z) to workbench orientation, given
-- array symmetry (any valid value of pa.symmetry).
--
-- The current coordinate system (coords) is one of
--
--  "mm"     -- workbench coordinates
--  "gu"     -- array volume coordinates
--  "abs_gu" -- absolute array coordinates (e.g removing symmetry)
--
-- (See section L.1.6 "Units and Coordinate Conventions" of
-- the printed SIMION manual for details.)
--
-- If "abs_gu", assumes (x,y,z) corresponds to a vector at point
-- (ion_px_gu,ion_py_gu,ion_pz_gu).
--
-- Note: mirroring value in symmetry parameter ignored.
--
-- This is designed to be called inside a SIMION accel_adjust segment.
local is_planar = {
  ["2dplanar"]=true,
  ["2dplanar[x]"]=true, ["2dplanar[y]"]=true, ["2dplanar[xy]"]=true,
  ["3dplanar"]="3d",
  ["3dplanar[x]"]="3d",  ["3dplanar[y]"]="3d",  ["3dplanar[z]"]="3d",
  ["3dplanar[xy]"]="3d", ["3dplanar[xz]"]="3d", ["3dplanar[yz]"]="3d",
  ["3dplanar[xyz]"]="3d"
}
local function to_wb_coords(symmetry, x,y,z, coords)
  -- Apply mirroring/symmetry to velocity to convert from PA coords
  -- to PA volume coords.
  -- This is the reverse of pa_coords_to_array_coords.
  if coords == 'abs_gu' then
    local planar = is_planar[symmetry]
    if planar then
      if ion_px_gu < 0 then x=-x end -- x mirrored
      if ion_py_gu < 0 then y=-y end -- y mirrored
      if planar == '3d' and ion_pz_gu < 0 then z=-z end -- z mirror (3D only)
    else -- cylindrical
      local r_gu = sqrt(ion_py_gu*ion_py_gu + ion_pz_gu*ion_pz_gu)
                 + 1.0e-10 -- (no div by zero)
      y, z = ion_py_gu / r_gu * y, ion_pz_gu / r_gu * y
      if ion_px_gu < 0 then x=-x end -- x mirrored
    end
    coords = 'gu'
  end
  if coords == 'gu' then
    -- Swing vector to workbench orientation.
    x,y,z = pa_orient_to_wb_orient(x,y,z)
  end
  return x,y,z
end


-- Updates local velocity with data at current ion position.
-- This is designed to be called from a SIMION accel_adjust or
-- other_actions segment.
-- Writes local_v[xyz]_mm_per_usec variables.
local function update_local_velocity(px,py,pz)
  -- Compute local velocity.
  local vel = vel_defs or velx_defs or vely_defs or velz_defs
  if vel then
    local vx,vy,vz
    if vel_defs then -- exists velocity vector function
      vx,vy,vz = vel_defs(px,py,pz)
    else             -- exists velocity vector components
      vx = velx_defs and velx_defs(px,py,pz) or 0
      vy = vely_defs and vely_defs(px,py,pz) or 0
      vz = velz_defs and velz_defs(px,py,pz) or 0
    end
    if vel_coords ~= 'mm' then
      vx,vy,vz = to_wb_coords(vel.symmetry, vx,vy,vz, vel_coords)  --rotate
    end
    local MM_USEC__M_SEC = 1.0e-3  -- (mm/usec)/(m/s)
    local f = SDS_v_field_scale_factor * MM_USEC__M_SEC
    local_vx_mm_per_usec = vx * f
    local_vy_mm_per_usec = vy * f
    local_vz_mm_per_usec = vz * f
  end
end


-- Updates local parameters (pressure, temperature, velocity, damping,
-- MFP, and average speed) with data at current ion position.
--
-- This is designed to be called from a SIMION accel_adjust segment.
--
-- Writes local_* variables and also variables written by update_ions_local.
--
local function update_local_parameters()
  -- Compute local potential gradient.
  -- This is used by both SDS_faims_mode == 1 and field-dependent mobility.
  if SDS_faims_mode ~= 0 or ions_alphabeta[ion_number] then
    -- Compute potential gradient (grad V) from acceleration (a) using
    -- Lorentz force law: m*a = F = q * -grad V.  Note that
    -- "a" (i.e. ion_a[xyz]_mm variables) actually incorporates all
    -- acceleration forces (charge repulsion, electric field and magnetic
    -- field), not just electric field (-ion_dvolts[xyz]_mm).  This
    -- acceleration is then re-expressed in terms of a potential gradient
    -- vector that can be fed into the mobility equation.
    -- Note that the charge repulsion acceleration is really just
    -- an electric field force.
    -- WARNING: acceleration vector here is non-relativistic.
    local M__MM = 1E-3
    local USEC2__SEC2 = 1E12
    local f = (ion_mass / ion_charge) * M__MM * M__MM * USEC2__SEC2
              * AMU_TO_KG / ELEMENTARY_CHARGE
    local_dvoltsx_Vmm = -ion_ax_mm * f
    local_dvoltsy_Vmm = -ion_ay_mm * f
    local_dvoltsz_Vmm = -ion_az_mm * f
    -- print('DEBUG-check', dvy_e[ion_number], ion_dvoltsy_mm)

    -- alternative:
    -- local_dvoltsx_Vmm = ion_dvoltsx_mm
    -- local_dvoltsy_Vmm = ion_dvoltsy_mm
    -- local_dvoltsz_Vmm = ion_dvoltsz_mm
  end

  -- Compute local pressure.
  if pres_defs then
    local px,py,pz
    if pres_coords == 'abs_gu' then
      px,py,pz = ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu
    elseif pres_coords == 'gu' then
      px,py,pz = ion_px_gu, ion_py_gu, ion_pz_gu
    else
      px,py,pz = ion_px_mm, ion_py_mm, ion_pz_mm
    end
    local_pressure_torr = pres_defs(px,py,pz) * abs(SDS_P_field_scale_factor)
    if local_pressure_torr == 0 then error("Array's local pressure = 0") end
  end

  -- Compute local temperature.
  if temp_defs then
    local px,py,pz
    if temp_coords == 'abs_gu' then
      px,py,pz = ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu
    elseif temp_coords == 'gu' then
      px,py,pz = ion_px_gu, ion_py_gu, ion_pz_gu
    else
      px,py,pz = ion_px_mm, ion_py_mm, ion_pz_mm
    end
    local_temperature_K = temp_defs(px,py,pz) * abs(SDS_T_field_scale_factor)
    if local_temperature_K == 0 then error("Array's local temperature = 0") end
  end

  -- Need to recalculate ion's local parameters.
  if pres_defs or temp_defs
     or ions_alphabeta[ion_number]  -- FAIMS
  then
    update_ions_local(ion_number)
  end

  -- Compute local temperature
  local px,py,pz
  if vel_coords == 'abs_gu' then
    px,py,pz = ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu
  elseif vel_coords == 'gu' then
    px,py,pz = ion_px_gu, ion_py_gu, ion_pz_gu
  else
    px,py,pz = ion_px_mm, ion_py_mm, ion_pz_mm
  end
  update_local_velocity(px,py,pz)
end


-- Apply (random-walk) diffusion effect to current ion's position
-- (ion_p*_mm) in this time-step delta_time (usec).
-- Given mass ratio ions_mr_ratio (mass_ion/mass_gas), mean-free-path
-- ions_MFP (mm), and mean thermal velocity ions_V (mm/usec) of ion.
-- Uses diffusion statistics (s_dist).
--
-- This is designed to be called inside a SIMION other_actions segment.
local function apply_diffusion(ions_mr_ratio, ions_MFP, ions_V, delta_time)
  -- Unitless distance traveled in N_DIST_COLLISIONS collisions,
  -- assuming normalized conditions.
  local dist_steps = diff_dist_steps(s_dist, ions_mr_ratio)

  -- Estimate number of collisions in current time step (delta_time),
  -- assuming local T and P (i.e. current MFP).
  local ncollisions =
      ions_V             -- (mm/usec) ion's current average thermal velocity
      / ions_MFP         -- (1/mm)    ion's current MFP
      * delta_time       -- (usec)    next time step

  -- Estimate jump distance (r) obtained over ncollisions collisions,
  -- assuming local T and P.  This is a simple scaling of dist_steps
  -- to actual conditions.
  local r =
    sqrt(ncollisions / s_dist.N_DIST_COLLISIONS) -- square root scaling law
    * dist_steps         -- r assuming MFP = 1.0 mm
    * ions_MFP           -- scale r with ion's MFP at local T and P

  -- Randomize direction of jump displacement vector (dx,dy,dz).
  local dx, dy, dz = rand()-0.5, rand()-0.5, rand()-0.5 -- (-0.5 to +0.5)
  local dr = sqrt(dx*dx + dy*dy + dz*dz)    -- original length of jump vector
  local f = r/dr                            -- scale to length r
  dx, dy, dz = dx*f, dy*f, dz*f

  -- Add jump vector to current ion position.
  ion_px_mm = ion_px_mm + dx
  ion_py_mm = ion_py_mm + dy
  ion_pz_mm = ion_pz_mm + dz
end


-- Apply Stokes' Law viscous damping (usec-1) in this time step
-- (ion_time_step) by damping acceleration (ion_a*_mm) given
-- mean bulk velocity of background gas (vx,vy,vz) (mm/usec).
--
-- This is designed to be called inside a SIMION accel_adjust segment.
local function apply_stokes_damping(damping, vx,vy,vz)
  -- See examples\drag\drag.lua for futher details on this implementation.
  if damping ~= 0 and ion_time_step ~= 0 then
    damping = abs(damping)  -- force positive

    local tterm = damping * ion_time_step  -- time constant
    local factor = (1 - exp(-tterm)) / tterm

    -- Store as new acceleration components.
    ion_ax_mm = factor*(ion_ax_mm - (ion_vx_mm - vx)*damping)
    ion_ay_mm = factor*(ion_ay_mm - (ion_vy_mm - vy)*damping)
    ion_az_mm = factor*(ion_az_mm - (ion_vz_mm - vz)*damping)
  end
end


-- BEGIN FAIMS
-- variables used by "quick RF cycle" feature
local last_time = {}
local last_x = {}
local last_y = {}
local last_z = {}
-- END FAIMS



-- SIMION segment called for each particle creation.
--
-- This is used to initialize some parameters.
local first = true
function M.segment.initialize()
  if SDS_enable == 0 then
    if first then first = false; print "SDS disabled" end
    return
  end

  -- Code executed on first call (basic init)
  if first then first = false
    if init then init() end

    if SDS_diffusion ~= 0 then
      print "SDS enabled (diffusion enabled)"
    else
      print "SDS enabled (diffusion disabled)"
    end
    -- Check adjustable variables.
    assert(SDS_pressure_torr ~= 0, "SDS_pressure_torr is zero")
    assert(SDS_temperature_K ~= 0, "SDS_temperature_K is zero")
    assert(SDS_collision_gas_mass_amu~=0,"SDS_collision_gas_mass_amu is zero")
    assert(SDS_collision_gas_diameter_nm ~= 0,
           "SDS_collision_gas_diameter_nm is zero")

    --seed(1) --DEBUG (disable randomization between runs)

    -- Set local temperature, pressure, and background gas velocity to
    -- global ones.  This is the default unless arrays/functions are
    -- defined.
    local_temperature_K = abs(SDS_temperature_K)
    local_pressure_torr = abs(SDS_pressure_torr)
    local MM_USEC__M_S = 1.0e-3  -- (mm/usec)/(m/s)
    local_vx_mm_per_usec = SDS_vx_m_per_sec * MM_USEC__M_S
    local_vy_mm_per_usec = SDS_vy_m_per_sec * MM_USEC__M_S
    local_vz_mm_per_usec = SDS_vz_m_per_sec * MM_USEC__M_S
    -- BEGIN FAIMS
    local_dvoltsx_Vmm = 0
    local_dvoltsy_Vmm = 0
    local_dvoltsz_Vmm = 0
    -- END FAIMS

    check_all_arrays()

    print_sds_parameters()

    -- Fill in mass data records, given adjusted vars on background gas.
    massdata_complete_records(massdata,
      SDS_collision_gas_diameter_nm, SDS_collision_gas_mass_amu)
  end -- first

  last_time = {} -- reset, FAIMS

  -- Get or estimate parameters for current ion.
  local mfpmass, vmass, ko,
        alpha, beta -- FAIMS
    = massdata_for_mass(massdata, ion_mass)

  -- Estimate Stokes' Law damping (usec-1)
  local emu = ELEMENTARY_CHARGE/AMU_TO_KG    -- (C kg-1 u)
  local damping =
      emu
      * 0.01       -- 10-4 * (10+6 usec sec-1)
      / ko         -- (10-4 m2 V-1 s-1) = 10-4 s C kg-1
      / ion_mass   -- (u)

  -- Compute log of mass ratio.
  local logmrratio = log10(ion_mass / SDS_collision_gas_mass_amu)

  -- Set current ion's parameters at standard temperature pressure (STP)
  -- or independent of position.
  local nion = ion_number
  ions_log_mr_ratio[nion]       = logmrratio
  -- BEGIN FAIMS
  ions_alpha[nion]              = alpha
  ions_beta[nion]               = beta
  ions_alphabeta[nion]          = alpha ~= 0 or beta ~= 0
  -- END FAIMS
  ions_STP_damping[nion]        = damping
  ions_STP_mfp_mm[nion]         = mfpmass
  ions_STP_Vo_mm_per_usec[nion] = vmass

  -- Compute initial values of current ion's parameters at current position.
  update_ions_local(ion_number)
end


-- BEGIN FAIMS
-- measured voltage gradient vector (V/mm) at end of time step (time_e usec).
-- these values are valid once all calls to accel_adjust in this time step
-- have completed.
local time_e = {}
-- END FAIMS


-- SIMION segment called to override time step sizes.
function M.segment.tstep_adjust()
  if SDS_enable == 0 then return end

  -- BEGIN FAIMS
  -- reset time_e before accel_adjust
  -- (note: tstep_adjust is called before accel_adjust).
  for k,v in next, time_e do time_e[k] = nil end
  -- END FAIMS

  -- Increase time step to minimum.
  ion_time_step = max(ion_time_step, abs(SDS_min_time_step_usec))
end


-- SIMION segment called to override acceleration vector.
--
-- This is used to apply the viscous mobility effect.
function M.segment.accel_adjust()
  if SDS_enable == 0 then return end

  -- WARNING: This segment can be limited to only certain PA
  -- instances.  It's particularly important that if magnetic and
  -- electrostatic arrays overlap then this segment should only
  -- be called in the magnetic one.
  if not instances[ion_instance] then return end

  -- Update local_* parameters.
  if is_local_params or SDS_faims_mode ~= 0
     or ions_alphabeta[ion_number] -- FAIMS
  then
    update_local_parameters()
  end

  -- BEGIN FAIMS
  -- Measure electric field vector at end of time-step.
  -- Note that accel_adjust is called over multiple points in the time
  -- step, but we only want the last point.
  -- The location in the time step can be particularly important when
  -- there is a sharp voltage transition (e.g. rectangular waveform edge)
  -- that coincides with the time step boundaries.
  if SDS_faims_mode ~= 0
     and ion_time_step >= (time_e[ion_number] or 0)
  then
    time_e[ion_number] = ion_time_step
    local_dvoltsx_last_Vmm[ion_number] = local_dvoltsx_Vmm
    local_dvoltsy_last_Vmm[ion_number] = local_dvoltsy_Vmm
    local_dvoltsz_last_Vmm[ion_number] = local_dvoltsz_Vmm
  end
  -- END FAIMS

  -- stokes' law viscous mobility effect.
  do
    if SDS_faims_mode ~= 0 then
      -- BEGIN FAIMS
      -- zero acceleration (terminal velocity)
      ion_ax_mm = 0; ion_ay_mm = 0; ion_az_mm = 0
      -- END FAIMS
    else
      local damping = ions_local_damping[ion_number]
      if not damping then
        error(string.format(
          "Particle #%d was not initialized.  Ensure that the particle " ..
          "originates strictly inside a potential array instance in which " ..
          "the initialize segment is called.", ion_number))
      end
      apply_stokes_damping(
        damping, local_vx_mm_per_usec,local_vy_mm_per_usec,local_vz_mm_per_usec)
    end
  end
end


-- SIMION segment called on every time-step.
--
-- This is used to apply the random-walk diffusion effect.
local first = true
function M.segment.other_actions()
  if SDS_enable == 0 then return end

  -- WARNING: This segment can be limited to only certain PA instances.
  -- It's particularly important that if magnetic and electrostatic arrays
  -- overlap, then this segment should only be called in the magnetic one.
  if not instances[ion_instance] then return end

  -- Code executed on first call.
  if first then first = false
    massdata_print(massdata, SDS_temperature_K, SDS_pressure_torr)
  end

  local delta_time = ion_time_step
  local allow_diffusion = true

  -- BEGIN FAIMS
  -- Force velocity of next time-step to be a direct function of force applied.
  if SDS_faims_mode ~= 0 then
    -- Note: this section uses the electric field at end of this time-step
    -- (rather than ion_dvoltsx_gu (yz), which represents the field at the
    -- beginning of the time step).  The field at the end of a voltage
    -- transition time-step reflects the field that exists in the next regular
    -- time step.

    -- Use local field and damping constant at end of time-step.
    local_dvoltsx_Vmm = local_dvoltsx_last_Vmm[ion_number]
    local_dvoltsy_Vmm = local_dvoltsy_last_Vmm[ion_number]
    local_dvoltsz_Vmm = local_dvoltsz_last_Vmm[ion_number]
    update_ions_local(ion_number)

    -- Compute ion mobility constant from local damping constant.
    local damping = ions_local_damping[ion_number]
    local emu = ELEMENTARY_CHARGE/AMU_TO_KG    -- (C kg-1 u)
    local ko = emu * 0.01 / damping / ion_mass -- (10-4 (m/s)/(V/m)) = 10-4 s C kg-1
    --print('DEBUG: local ko=', ko)

    -- Set velocity vector (relative to backgroung gas bulk velocity) for next
    -- time-step according to mobility equation v = ko * E.
    local ko_conversion = ko * (1000 * 1E-4 * 1E-3)
    ion_vx_mm = local_vx_mm_per_usec - local_dvoltsx_Vmm * ko_conversion
    ion_vy_mm = local_vy_mm_per_usec - local_dvoltsy_Vmm * ko_conversion
    ion_vz_mm = local_vz_mm_per_usec - local_dvoltsz_Vmm * ko_conversion
  end

  -- Optional "quick RF cycle" optimization.  See comments
  -- above SDS_quick_cycle adjustable variable for description.
  if SDS_quick_cycles ~= 0 then
    -- Begin normal RF cycle, if not yet started,
    -- and record conditions at start of cycle.
    if not last_time[ion_number] then
      last_time[ion_number] = ion_time_of_flight
      last_x[ion_number] = ion_px_mm
      last_y[ion_number] = ion_py_mm
      last_z[ion_number] = ion_pz_mm
    end

    -- Disable diffusion during normal RF cycle trace.
    -- (Diffusion can be applied later.)
    -- Note that displacements calculated in the normal RF cycle are repeated in
    -- each of the quick RF cycles.  If we enabled randomized diffusion effects
    -- in the normal RF cycle, those effects in the quick RF cycles will not be
    -- randomized and independent.  Diffusion distance would scale linearly with
    -- time rather than by the "square root scaling law" (random walk)
    -- discussed elsewhere.
    allow_diffusion = false

    -- Detect end of normal RF cycle.
    local dt = ion_time_of_flight - last_time[ion_number]
    local dtp = dt / SDS_quick_period
    local dtq = abs(dtp - 1)
    if dtq < 0.001 then -- i.e. very close
      -- compute displacement over RF cycle.
      local dx = ion_px_mm - last_x[ion_number]
      local dy = ion_py_mm - last_y[ion_number]
      local dz = ion_pz_mm - last_z[ion_number]

      local NQUICK
      if SDS_quick_cycles == -1 then -- automatically choose
        -- Set number of cycles to as to advance approximately 1 grid unit.
        local dru = sqrt(dx*dx+dy*dy+dz*dz) / ion_mm_per_grid_unit
        local MAX_QUICK_CYCLES = 1000
        NQUICK = min(floor(1 / dru), MAX_QUICK_CYCLES)
        --print('DEBUG:nquick', NQUICK)
      else
        NQUICK = SDS_quick_cycles
      end

      -- Simulate quick RF cycles by jumping ion.
      ion_time_of_flight = ion_time_of_flight + dt * NQUICK
      ion_px_mm = ion_px_mm + dx * NQUICK
      ion_py_mm = ion_py_mm + dy * NQUICK
      ion_pz_mm = ion_pz_mm + dz * NQUICK
      last_time[ion_number] = nil  -- clear normal RF cycle

      -- Now apply diffusion over delta_time (normal cycle + quick RF cycles)
      allow_diffusion = true
      delta_time = dt * (NQUICK+1)
    elseif dtp > 1 then -- failed to accurately measure normal RF cycle
                       -- (e.g. time steps not accurately aligned to RF period)
      last_time[ion_number] = nil  -- safely skip quick RF cycle optimization
      allow_diffusion = true -- don't forget diffusion
      delta_time = dt
    end
  end

  -- END FAIMS

  -- random-walk diffusion effect.
  if SDS_diffusion ~= 0 and allow_diffusion then
    apply_diffusion(ions_log_mr_ratio[ion_number],
                    ions_local_mfp_mm[ion_number],
                    ions_local_V_mm_per_usec[ion_number],
                    delta_time)
  end
  --print('time constant (usec)=', 1 / ions_local_damping[ion_number])
  -- sim_update_pe_surface = 1
  -- print_coordinates()
end


local function merge_segments(t)
  for name,newseg in pairs(t) do
    local oldseg = segment[name]
    segment[name] =
      oldseg and function() oldseg(); newseg() end
             or  newseg
  end
end


-- Install SIMION segments (copy functions from M.segment table into
-- segment table).
function M.install()
  merge_segments(M.segment)
end


-- Install segments on load.
if opt ~= 'noinstall' then
  M.install()
end


-- Helper function.
local function guess_vel_coords(v)
  -- If v is an array (or SIMION PA), assume velocity vectors in array
  -- are oriented with respect to the array (absolute array
  -- coordinates).  Otherwise (e.g. it is a
  -- function), assume oriented with respect to the workbench
  -- (workbench coordinates).
  local is_array = type(v) == 'table' or type(v) == 'userdata' and v.symmetry
  vel_coords = is_array and 'abs_gu' or 'mm'
end
local function guess_pres_coords(v)
  local is_array = type(v) == 'table' or type(v) == 'userdata' and v.symmetry
  pres_coords = is_array and 'abs_gu' or 'mm'
end
local function guess_temp_coords(v)
  local is_array = type(v) == 'table' or type(v) == 'userdata' and v.symmetry
  temp_coords = is_array and 'abs_gu' or 'mm'
end


-- Catch getting and setting on this module table.
-- For example, support things like
--   local SDS = simion.import("collision_sds.lua")
--   SDS.pressure = function() return ion_px_mm * 2 end
local mt = {}
function mt.__index(t,k)  -- GET t[k]
  if     k == 'pressure'    then return pres_defs
  elseif k == 'temperature' then return temp_defs
  elseif k == 'velocity'    then return vel_defs
  elseif k == 'velocity_x'  then return velx_defs
  elseif k == 'velocity_y'  then return vely_defs
  elseif k == 'velocity_z'  then return velz_defs
  elseif k == 'velocity_coordinates'
                            then return vel_coords
  elseif k == 'instances'   then return instances
  elseif k == 'init'        then return init
  else error(k) end
end
function mt.__newindex(t,k,v)  -- SET t[k]=v
  print("SDS Defining " .. k .. ":", v)
  if     k == 'pressure'    then pres_defs = v; guess_pres_coords(v)
  elseif k == 'temperature' then temp_defs = v; guess_temp_coords(v)
  elseif k == 'velocity'    then vel_defs  = v; guess_vel_coords(v)
  elseif k == 'velocity_x'  then velx_defs = v; guess_vel_coords(v)
  elseif k == 'velocity_y'  then vely_defs = v; guess_vel_coords(v)
  elseif k == 'velocity_z'  then velz_defs = v; guess_vel_coords(v)
  elseif k == 'velocity_coordinates'
                            then vel_coords = v
  elseif k == 'instances'   then instances = v
  elseif k == 'init'        then init = v
  else error(k) end
end
setmetatable(M, mt)


return M
