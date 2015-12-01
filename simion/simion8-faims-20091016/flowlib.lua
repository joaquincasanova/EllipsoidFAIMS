-- flowlib.lua
-- Gas flow modeling library.
-- Calculates viscosity constants and Poiseuille flows.
--
-- D.Manura. 2009
-- (c) 2009 Scientific Instrument Services, Inc. (Licensed
-- under SIMION 8.1).


-- detect program errors ( http://simion.com/issue/490 )
if checkglobals then checkglobals() end

local M = {}


-- Utility functions to check function parameter tables.
-- You can mostly ignore this.
local _T = {}
local _CHECKED = {}
local function begin_check(t)
  local checked = {}
  local mt = {}
  function mt.__index(tw, k)
    checked[k] = true
    local v = t[k]
    if v == nil then
      error('field "' .. k .. '" not found')
    end
    return v
  end
  local tw = {[_T] = t, [_CHECKED] = checked}
  return setmetatable(tw, mt)
end
local function end_check(tw)
  local t = tw[_T]
  local checked = tw[_CHECKED]
  for k in pairs(t) do
    if not checked[k] then
      error('field "' .. k .. '" not expected')
    end
  end
  return t
end


-- Gas properties used to compute viscosity of
-- a gas at given temperature via Sutherland's formula.
--
-- C = Sutherland's constant (K)
-- T0 = reference temperature (degrees K)
-- mu0 = reference viscosity (Pa s)
--
-- From http://www.lmnoeng.com/Flow/GasViscosity.htm
-- based on Crane (1988, p.A-5) and CRC (1984, pp.F-42-44).
local gas_properties = {
  air = {C = 120, T0 = 524.07*(5/9), mu0 = 0.01827E-3},
  NH3 = {C = 370, T0 = 527.67*(5/9), mu0 = 0.00982E-3},
  CO2 = {C = 240, T0 = 527.67*(5/9), mu0 = 0.01480E-3},
  CO  = {C = 118, T0 = 518.67*(5/9), mu0 = 0.01720E-3},
  H2  = {C = 72,  T0 = 528.93*(5/9), mu0 = 0.00876E-3},
  N2  = {C = 111, T0 = 540.99*(5/9), mu0 = 0.01781E-3},
  O2  = {C = 127, T0 = 526.05*(5/9), mu0 = 0.02018E-3},
  SO2 = {C = 416, T0 = 528.57*(5/9), mu0 = 0.01254E-3}
}


-- Compute mobility of gas at given temperature.
-- gas is string containing gas name in gas_properties
-- (e.g. "air") or a table in the format of gas_properties["air"].
-- T is temperature in Kelvin (K).
-- return mobility constant in Pa s.
-- Uses Sutherland's formula:
--   mu = mu0 * ((T0+C)/(T+C)) * (T/T0)^(3/2)
function M.compute_mu(gas, T)
  local props = type(gas) == 'string' and gas_properties[gas] or gas
  if not props then
    error("properties for gas " .. gas .. " not found")
  end
  local C   = props.C
  local T0  = props.T0
  local mu0 = props.mu0
  local mu = mu0 * ((T0 + C)/(T + C)) * (T/T0)^(3/2)
  return mu
end


-- Defines Poiseuille [4] gas flow between infinite parallel plates
-- inside SDS.
-- The velocity is parabolic, with maximum velocity (vx_m_psec)
-- in the gap center and zero velocity at the plates.
-- Pressure is assumed to be p0_torr at x=0 mm
-- and decreases for increasing x according to the
-- Poiseuille planar flow equation.
--
-- t is a table containing these fields:
--
--   SDS  - SDS object.
--   d_mm - distance between plates (mm)
--   x0_mm     - x origin (mm) where p0 is measured
--   y0_mm     - y center (mm) half-way between plates
--   p0_torr   - pressure at x=x0_mm (torr)
--   vx_m_psec - maximum gas velocity, in center of plates.
--               This is (3/2) the average velocity.
--   mu_pa_s -  dynamic viscosity (Pa s).
--              Note: approx 1.861E-5 for air at 298 K.
--
function M.define_poiseuille_planar(t)
  t = begin_check(t)
  local SDS          = t.SDS
  local d_mm         = t.d_mm
  local x0_mm        = t.x0_mm
  local y0_mm        = t.y0_mm
  local p0_torr      = t.p0_torr
  local vx_m_psec    = t.vx_m_psec
  local mu_pa_s      = t.mu_pa_s
  t = end_check(t)

  local MM_M    = 1000       -- mm/m
  local TORR_PA = 760/101325 -- Torr/Pa

  SDS.velocity = function(x,y,z)
    -- Parabolic profile.
    local r = (y - y0_mm) / (0.5*d_mm)  -- (unitless, normalized to 1)
    local vx = vx_m_psec * (1 - r*r)    -- (m/sec)
    return vx, 0, 0
  end

  SDS.pressure = function(x,y,z)
    local dp_dx_torr_mm =     -- pressure gradient (torr/mm)
      vx_m_psec * (-8) * mu_pa_s * (d_mm/MM_M)^-2   -- Pa/m
       * TORR_PA / MM_M
    local p = p0_torr + (x-x0_mm) * dp_dx_torr_mm
    --print('DEBUG',p, dp_dx_torr_mm)
    assert(p >= 0, "pressure became negative")
    return p
  end

  -- note: SDS.temperature not customized here.
end


-- Defines Poiseuille [4] gas flow between infinite concentric cylinders
-- inside SDS.  Approximates gap as infinite parallel plates.
--
-- Note: similar to planar version in that
-- concentric cylinder plates locally resemble parrallel plates.
--
-- t is a table containing these fields:
--
--  SDS  - SDS object.
--  d_mm      -- Distance between plates (mm)
--  x0_mm     -- X origin (mm)
--  y0_mm     -- Y center of cylinders (mm)
--  z0_mm     -- Z center of cylinders (mm)
--  p0_torr   -- Pressure at x=x0_mm (torr)
--  vx_m_psec -- maximum gas velocity, in center of plates.
--               This is (3/2) the average velocity.
--  mu_pa_s   -- dynamic viscosity (Pa s)
--
function M.define_poiseuille_cylindrical(t)
  t = begin_check(t)
  local SDS       = t.SDS
  local d_mm      = t.d_mm
  local x0_mm     = t.x0_mm
  local y0_mm     = t.y0_mm
  local z0_mm     = t.z0_mm
  local p0_torr   = t.p0_torr
  local vx_m_psec = t.vx_m_psec
  local mu_pa_s   = t.mu_pa_s
  t = end_check(t)

  local MM_M = 1000  -- mm/m
  local TORR_PA = 760/101325 -- Torr/Pa

  SDS.velocity = function(x,y,z)
    -- Parabolic profile.
    local R0 = 1.25  -- radius of mipoint between cylinders (mm)
    local r = (sqrt((y-y0_mm)^2 + (z-z0_mm)^2) - R0)
              / (0.5*d_mm)  -- (unitless, normalized to 1)
    local vx = vx_m_psec  -- (m/s)
    return vx * (1 - r*r), 0, 0
  end

  SDS.pressure = function(x,y,z)
    local dp_dx_torr_mm =     -- pressure gradient (torr/mm)
       vx_m_psec * (-8) * mu_pa_s * (d_mm/MM_M)^-2   -- Pa/m
       * TORR_PA / MM_M
    local p = p0_torr + x * dp_dx_torr_mm
    --print('DEBUG',p, dp_dx_torr_mm)
    assert(p >= 0, "pressure became negative")
    return p
  end

  -- note: SDS.temperature not customized here.
end


return M
