-- spectrum_waveform.lua - SIMION workbench user program
-- Aquires a FAIMS spectrum .  The spectrum (intensity v.s. CV) data
-- is sent to the Log and also plotted in Excel.  Also uses the
-- waveform library (waveformlib.lua).
--
-- This differs from spectrum.lua in that it defines the waveform
-- using the waveform library (waveformlib.lua).  The results should
-- be the same.
--
-- WARNING: waveformlib.lua currently does not support SDS_faims_mode.
-- 
-- D.Manura, 2008-07.

simion.workbench_program()

-- Load SDS collisional model.
local SDS = simion.import 'collision_sds.lua'

-- Load waveform library
local WAVE = simion.import 'waveformlib.lua'

-- Load Excel plotting library.
local EXCEL = simion.import 'excellib.lua'


--## BEGIN USER ADJUSTABLE VARIABLES

  -- -- 1 = enable FAIMS optimizations (normally don't change this)
  -- adjustable SDS_faims_mode = 1

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

--## END USER ADJUSTABLE VARIABLES


-- Compensation voltage
-- This is controlled by the program.
local wave_CV = 0

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

-- Install waveform library segments.
WAVE.install {
  -- Update PE surface display every this number of microseconds.
  -- This is optional.  Remove or set to nil to disable.
  pe_update_period = 1;
}

-- SIMION segment called upon each particle initialization.
local old_initialize = segment.initialize
function segment.initialize()
  old_initialize()

  nions = ion_number  -- note: last write contains last ion number

  if ion_number == 1 then  -- first ion
    nrun = nrun + 1

    -- Set up parameters for this run.
    wave_CV = SPECTRUM_CV_begin + (nrun-1) * SPECTRUM_CV_step
    intensity = 0   -- reset

    -- Set up waveform.
    ---- low voltage in cycle.
    local Vlow = - wave_DV * (wave_duty / (1 - wave_duty))
    WAVE.set_waveform(
      -- Define waveform for each adjustable electrode.
      -- Note: times in microseconds.
      WAVE.waveforms {
        WAVE.electrode(2) {
          WAVE.loop(math.huge) {
            WAVE.lines {
              {time=0,                       potential=wave_CV + wave_DV};
              {time=wave_period/256,         potential=wave_CV + wave_DV};
              {time=wave_period * wave_duty, potential=wave_CV + wave_DV};
              {time=wave_period * wave_duty, potential=wave_CV + Vlow};
              {time=wave_period * (wave_duty+1/256), potential=wave_CV + Vlow};
              {time=wave_period,             potential=wave_CV + Vlow};
            };
          };
        };
      }
    )
  end
end


-- SIMION segment called on each time-step for each particle.
local old_other_actions = segment.other_actions
function segment.other_actions()
  old_other_actions()

  -- Detect ions (you may want to modify how this works).
  if ion_px_mm > x_detect then
    intensity = intensity + 1
    ion_splat = 1  -- splat it
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
    print("run=,", nrun, ",CV=,", wave_CV, ",intensity=,", intensity)
    spectrum[#spectrum+1] = {wave_CV, intensity}

    if sim_rerun_flym == 0 then  -- last run
      EXCEL.plot(spectrum)   -- plot spectrum in Excel.
    end
  end
end

local old_tstep_adjust = segment.tstep_adjust
function segment.tstep_adjust()
  old_tstep_adjust()

  -- fraction of current pulse.
  local f = (ion_time_of_flight / wave_period) % 1

--FIX
  if f == 0 or f == wave_duty then
     -- mark()
    _G.transition = true
  else
    _G.transition = false
  end
end
