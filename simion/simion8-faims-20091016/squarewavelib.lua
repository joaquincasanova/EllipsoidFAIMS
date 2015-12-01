-- squarewavelib.lua
-- Simple square waveform generation.

local M = {}

-- Period of waveform (usec)
adjustable wave_period = 2.0
-- Fraction of period waveform is at high voltage.
adjustable wave_duty = 0.25
-- Dispersion voltage
adjustable wave_DV = 750
-- Compensation voltage
adjustable wave_CV = 0
-- Minimum number of time-steps per waveform period.
-- The program may use a larger value.
-- Set to 0 to have this controlled instead by SIMION via
-- the trajectory quality (TQual) parameter.
adjustable wave_timesteps = 1

-- From SDS model.
adjustable SDS_quick_period

function M.fast_adjust()
  local Vlow = - wave_DV * (wave_duty / (1 - wave_duty))
  local f = ((ion_time_of_flight) / wave_period) % 1
  local V = f < wave_duty and wave_DV or Vlow
  --print('V=',V, 't=', ion_time_of_flight, ion_time_step, ion_px_mm)
  adj_elect[2] = wave_CV + V
end

-- Improves time-step accuracy, causing time-steps to to end closer to
-- pulse edges.
local min = math.min
function M.tstep_adjust()
  -- In case quick RF cycle feature used, ensure period is properly set.
  SDS_quick_period = wave_period

  -- initially, increase time-step to a fraction of the RF period.
  if wave_timesteps ~= 0 then
    ion_time_step = wave_period / wave_timesteps
  end

  -- prevents lockups due to numerical error when
  -- 0 < ion_time_step < TOL and
  -- ion_time_of_flight + ion_time_step == ion_time_of_flight.
  -- In units of microseconds.
  local TOL = 1e-10

  -- A very small time-step will cover the region of the voltage transition by
  -- this fraction of wave_period before and after the voltage transition.
  -- Isolating voltage transitions in this small time step allows the regular
  -- time steps to have constant voltages, thereby avoiding possible errors
  -- and complications changing voltages during time-steps.
  local TOL2 = 1e-6

  -- voltage transition points as a fraction of wave_period.
  local f1 = TOL2
  local f2 = wave_duty - TOL2
  local f3 = wave_duty + TOL2
  local f4 = 1 - TOL2

  -- current fraction of wave period.
  local f = (ion_time_of_flight / wave_period) % 1

  -- Reduce time step to hit next transition point.
  -- Set the flag _G.transition to true iff this is
  -- one of the tiny voltage transition time steps.
  if f < f1 - TOL then
    local next = (f1 - f) * wave_period
    ion_time_step = min(ion_time_step, next)
  elseif f < f2 - TOL then
    local next = (f2 - f) * wave_period
    ion_time_step = min(ion_time_step, next)
  elseif f < f3 - TOL then
    local next = (f3 - f) * wave_period
    ion_time_step = min(ion_time_step, next)
  elseif f < f4 - TOL then
    local next = (f4 - f) * wave_period
    ion_time_step = min(ion_time_step, next)
  else
    local next = (1 + f1 - f) * wave_period
    ion_time_step = min(ion_time_step, next)
  end

  -- mark()
end


-- install
local old = segment.fast_adjust
function segment.fast_adjust()
  if old then old() end
  M.fast_adjust()
end
local old = segment.tstep_adjust
function segment.tstep_adjust()
  if old then old() end
  M.tstep_adjust()
end
