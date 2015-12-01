-- cylindrical_plates.lua - SIMION workbench user program.
-- Demonstrates FAIMS using an ideal cylindrical FAIMS cell.
--
-- D.Manura, 2008-06

simion.workbench_program()

-- Load waveform library (import only one)
simion.import 'squarewavelib.lua'
--simion.import 'bisinusoidalwavelib.lua'

-- Load SDS collisional model.
local SDS = simion.import 'collision_sds.lua'


-- 1 = enable FAIMS optimizations (normally don't change this)
adjustable SDS_faims_mode = 1

-- Period of waveform (usec)
adjustable wave_period = 2.0
-- Fraction of period waveform is at high voltage.
adjustable wave_duty = 0.25
-- Dispersion voltage
adjustable wave_DV = 750
-- Compensation voltage
adjustable wave_CV = 0
-- Minimum number of time-steps per waveform period.
-- (note: larger value required due to non-uniform fields in cylindrical gap)
adjustable wave_timesteps = 256

-- (m/s)
adjustable SDS_vx_m_per_sec = 6.667

-- Whether to enable SDS randomized diffusion effect.  1=enabled,
-- 0=disabled.  WARNING! Normally, you want this enabled (1) for the
-- most realistic calculation.  However, for comprehension purposes,
-- it may be useful to temporarily disable (0) it so that only the
-- SDS mobility effect (which facilitates FAIMS) is present.
adjustable SDS_diffusion = 0

-- Terminate particles early if particle x position (mm)
-- is greater than this value.  (Speeds simulation.)
adjustable x_max = 3


-- SIMION segment called on each time step.
local old = segment.other_actions
function segment.other_actions()
  old()
  -- mark()

  if ion_px_mm > x_max then
  -- if ion_time_of_flight > 133333.3333333 then
    print('Note: User program terminating particle early to speed simulation.')
    ion_splat = -4
  end
end


if SDS_diffusion == 0 then
  print 'Warning: diffusion effect is by default disabled.'
  print '  Set SDS_diffusion to 1 to enable.'
end


-- Poiseuille planar flow between concentric cylinders.
-- Remove following "--[[" to enable this code.

function SDS.init()
  local FLOW = simion.import "flowlib.lua"
  adjustable SDS_pressure_torr
  adjustable SDS_temperature_K
  adjustable SDS_collision_gas_mass_amu
  assert(SDS_collision_gas_mass_amu == 28.94515,
    "FLOW assumes air but SDS_collision_gas_mass_amu not 28.94515)")
  local mu_pa_s = FLOW.compute_mu('air', SDS_temperature_K)
  print('mu=', mu_pa_s, 'Pa s')
  FLOW.define_poiseuille_cylindrical {
    SDS = SDS,
    d_mm = 0.5,  -- Distance between plates (mm)
    x0_mm = 0,      -- X origin (mm)
    y0_mm = -1.25,  -- Y center of cylinders (mm)
    z0_mm = 0,      -- Z center of cylinders (mm)
    p0_torr = SDS_pressure_torr,  -- Pressure at origin (torr)
    vx_m_psec = SDS_vx_m_per_sec, -- max velocity (m/s) in center
    mu_pa_s = mu_pa_s  -- dynamic viscosity (Pa s)
  }
end

