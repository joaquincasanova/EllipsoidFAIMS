-- plates_waveform.lua - SIMION workbench user program.
-- Demonstrates FAIMS using an ideal planar FAIMS cell,
-- using the waveform library (waveformlib.lua).
--
-- This differs from plates.lua in that it defines the waveform using
-- the waveform library (waveformlib.lua).  The results are expected
-- to be the same as in plates.lua.
--
-- WARNING: waveformlib.lua currently does not support SDS_faims_mode.
--
-- D.Manura, 2008-07

simion.workbench_program()

-- Load waveform library.
local WAVE = simion.import 'waveformlib.lua'

-- Load SDS collisional model.
local SDS = simion.import 'collision_sds.lua'


--## BEGIN USER ADJUSTABLE VARIABLES

  -- WARNING! NOT ENABLED
  -- -- 1 = enable FAIMS optimizations (normally don't change this)
  -- adjustable SDS_faims_mode = 1

  -- Period of waveform (usec)
  adjustable wave_period = 2.0
  -- Fraction of period waveform is at high voltage.
  adjustable wave_duty = 0.25
  -- Dispersion voltage
  adjustable wave_DV = 750

  -- 1 L/min through cross-sectional area of 5 * 0.5 mm^2.
  adjustable SDS_vx_m_per_sec = 6.667

  -- Whether to enable SDS randomized diffusion effect.  1=enabled,
  -- 0=disabled.  WARNING! Normally, you want this enabled (1) for the
  -- most realistic calculation.  However, for comprehension purposes,
  -- it may be useful to temporarily disable (0) it so that only the
  -- SDS mobility effect (which facilitates FAIMS) is present.
  adjustable SDS_diffusion = 0

  -- Whether to plot waveform in Excel (1=yes,0=no).
  adjustable plot_enable = 0

  -- Maximum time value in plot (usec) (math.huge=+infinity)
  adjustable plot_maxtime = 2

--## END USER ADJUSTABLE VARIABLES


-- Compensation voltage
-- This is controlled by the program.
local wave_CV = 0


-- Install waveform library segments.
WAVE.install {
  -- Update PE surface display every this number of microseconds.
  -- This is optional.  Remove or set to nil to disable.
  pe_update_period = 1;
}


-- SIMION segment called on each particle initialization in PA.
local old_initialize = segment.initialize
function segment.initialize()
  old_initialize()
  if ion_number == 1 then
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
              {time=wave_period * wave_duty, potential=wave_CV + wave_DV};
              {time=wave_period * wave_duty, potential=wave_CV + Vlow};
              {time=wave_period,             potential=wave_CV + Vlow};
            };
          };
        };
      }
    )

    -- Plot waveform (this is optional and can be removed).
    if plot_enable ~= 0 then WAVE.plot_waveform(nil, plot_maxtime) end
  end
end


-- SIMION segment called on each time step.
-- WARNING!
--   Terminates particles prematurely to speed up the simulation.
--   You could remove this.
local old_other_actions = segment.other_actions
function segment.other_actions()
  old_other_actions()
  if ion_px_mm > 1 then ion_splat = 1 end
end


