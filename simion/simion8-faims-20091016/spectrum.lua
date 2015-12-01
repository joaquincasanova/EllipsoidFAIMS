-- spectrum.lua - SIMION workbench user program
-- Aquires a FAIMS spectrum.  The spectrum (intensity v.s. CV) data
-- is sent to the Log and also plotted in Excel.
--
-- D.Manura, 2008-07.

simion.workbench_program()


-- Load waveform library.
local WAVE = simion.import 'squarewavelib.lua'

-- Load SDS collisional model.
local SDS = simion.import 'collision_sds.lua'


-- Load Excel plotting library.
local EXCEL = simion.import 'excellib.lua'


--## BEGIN USER ADJUSTABLE VARIABLES

  -- 1 = enable FAIMS optimizations (normally don't change this)
  adjustable SDS_faims_mode = 1

  -- 1 L/min through cross-sectional area of 5 * 0.5 mm^2.
  adjustable SDS_vx_m_per_sec = 6.667

  -- 0 = disable diffusion (makes FAIMS pattern clearer)
  -- 1 = enable diffusion  (normal use)
  adjustable SDS_diffusion = 1

  -- FAIMS spectrum scan defined from begin to end voltages in steps of
  -- step voltage.
  adjustable SPECTRUM_CV_begin = -10
  adjustable SPECTRUM_CV_end = 10
  adjustable SPECTRUM_CV_step = 0.25

  -- Particles that reach this x (mm) are considered detected
  -- by the detector.  You may want to modify how this works.
  adjustable x_detect = 5

  -- Period of waveform (usec)
  adjustable wave_period = 2.0
  -- Fraction of period waveform is at high voltage.
  adjustable wave_duty = 0.25
  -- Dispersion voltage
  adjustable wave_DV = 750

  -- Whether to plot results in Excel (0=no, 1=yes)
  adjustable excel_enable = 1

  adjustable workaround_max_y = 0.259

--## END USER ADJUSTABLE VARIABLES
if SDS_diffusion == 0 then
  print('Warning: diffusion effect is by default disabled.  Set SDS_diffusion to 1 to enable.')
end


-- Compensation voltage
-- This is controlled by the program.
adjustable wave_CV = 0

-- Current run number.
local nrun = 0

-- Particle count (scan intensity) for current scan.
local intensity

-- Total number of ions.
local nions

-- Spectrum.  This is stored as a table plottable with the Excel
-- library.
local spectrum = {}
spectrum.header = {'CV', 'intensity'}
spectrum.title = 'FAIMS Spectrum'
spectrum.lines = true


-- SIMION segment called upon each particle initialization.
local old = segment.initialize
function segment.initialize()
  old()

  nions = ion_number  -- note: last write contains last ion number

  if ion_number == 1 then  -- first ion
    nrun = nrun + 1

    -- Set up parameters for this run.
    wave_CV = SPECTRUM_CV_begin + (nrun-1) * SPECTRUM_CV_step
    intensity = 0   -- reset
  end
end


-- SIMION segment called on each time-step for each particle.
local old = segment.other_actions
function segment.other_actions()
  old()

  -- Detect ions (you may want to modify how this works).
  if ion_px_mm > x_detect then
    intensity = intensity + 1
    ion_splat = 4  -- splat it
  end

  -- Workaround problem where if SDS jumps particle outside array,
  -- terminate segment is not called. (TO DO: improve SIMION to avoid this)
  if ion_py_mm >= workaround_max_y then
    ion_py_mm = workaround_max_y
    ion_splat = 4
  end
  if ion_py_mm <= -workaround_max_y then
    ion_py_mm = -workaround_max_y
    ion_splat = 4
  end
end


-- SIMION segment called for each particle following termination
-- of all particles.
assert(not segment.terminate)
function segment.terminate()

  --FIX: warning: if particles splat outside the PA instance, this
  --doesn't get called.

  if ion_number == nions then  -- last ion
    -- Determine whether to do another run.
    sim_rerun_flym =
      (wave_CV + SPECTRUM_CV_step <= SPECTRUM_CV_end) and 1 or 0

    -- Record and store current scan result in spectrum.
    print("run=", nrun, " CV=", wave_CV, ",intensity=,", intensity)
    spectrum[#spectrum+1] = {wave_CV, intensity}

    if sim_rerun_flym == 0 then  -- last run
      if excel_enable ~= 0 then
        EXCEL.plot(spectrum)   -- plot spectrum in Excel.
      end
    end
  end
end


-- Poiseuille planar flow.
-- Remove following "--[[" to enable this code.
--[[
function SDS.init()
  local FLOW = simion.import "flowlib.lua"
  adjustable SDS_pressure_torr
  adjustable SDS_temperature_K
  adjustable SDS_collision_gas_mass_amu
  assert(SDS_collision_gas_mass_amu == 28.94515,
    "FLOW assumes air but SDS_collision_gas_mass_amu not 28.94515.")
  local mu_pa_s = FLOW.compute_mu('air', SDS_temperature_K)
  print('mu=', mu_pa_s, 'Pa s')
  FLOW.define_poiseuille_planar {
    SDS = SDS,
    d_mm = 0.5,   -- Distance between plates (mm)
    x0_mm = 0,    -- X origin (mm)
    y0_mm = 0,    -- Y center (mm)
    p0_torr = SDS_pressure_torr,  -- Pressure at origin (torr)
    vx_m_psec = SDS_vx_m_per_sec, -- max velocity (m/s) in center
    mu_pa_s = mu_pa_s  -- dynamic viscosity (Pa s)
  }
end
--]]
