-- bisinusoidalwavelib.lua
-- Simple bisinusoidal waveform generation.

local M = {}

-- Period of waveform (usec)
adjustable wave_period = 2.0
-- Dispersion voltage
adjustable wave_DV = 750
-- Compensation voltage
adjustable wave_CV = 0
-- Minimum number of time-steps per waveform period.
-- The program may use a larger value.
-- Set to 0 to have this controlled instead by SIMION via
-- the trajectory quality (TQual) parameter.
adjustable wave_timesteps = 16

-- From SDS model.
adjustable SDS_quick_period

-- constants for bisinusoidal shape
local h = 2
local F = 2

local PI = math.pi

function M.fast_adjust()
  local W = (1/wave_period) * 2 * PI
  local t = ion_time_of_flight
  adj_elect[2] = wave_CV + (F*sin(W*t) + sin(h*W*t - 0.5*PI))*wave_DV/(F + 1)

  -- print('V=',adj_elect[4], 't=', ion_time_of_flight, ion_time_step,ion_px_mm)
end

-- Time-step control.  Must be a sufficiently small fraction
-- of period to fully represent wave form.
local min = math.min
function M.tstep_adjust()
  -- In case quick RF cycle feature used, ensure period is properly set.
  SDS_quick_period = wave_period

  ion_time_step = wave_period / wave_timesteps
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
