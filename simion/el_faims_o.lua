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

-- Load waveform library.
local WAVE = simion.import 'squarewavelib.lua'

-- Load Excel plotting library.
local EXCEL = simion.import 'excellib.lua'

filename = "el_faims.csv"

filehandle = io.output(filename)

postot = 0

--## BEGIN USER ADJUSTABLE VARIABLES

  -- -- 1 = enable FAIMS optimizations (normally don't change this)
  adjustable SDS_faims_mode = 1

  --.
  --adjustable SDS_vz_m_per_sec = .15

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
  adjustable z_detect = 14.7
  adjustable x_detect = 0
  adjustable y_detect = 0
  adjustable r_detect = 1
   
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
local old_other_actions = segment.other_actions
function segment.other_actions()
  old_other_actions()

  -- Detect ions (you may want to modify how this works).
  if ion_pz_mm > z_detect then
	local pos = (ion_px_mm-x_detect)^2+(ion_py_mm-y_detect)^2
	--print(pos)
	if pos <=r_detect^2 then
		intensity = intensity + 1
		postot = sqrt(pos)+postot
		ion_splat = 1  -- splat it
	end
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
 
    print(nrun, wave_CV, intensity,postot/intensity)
	io.write(nrun, ", ", wave_CV,", ", intensity,", ",postot/intensity,"\n")
  end
  
  if sim_rerun_flym == 0 then  -- last run
	io.flush()
  end
end

-- Poiseuille planar flow.
-- Remove following "--[[" to enable this code.

function SDS.init()
  local FLOW = simion.import "flowlib.lua"
  adjustable SDS_pressure_torr
  adjustable SDS_temperature_K
  adjustable SDS_collision_gas_mass_amu
  assert(SDS_collision_gas_mass_amu == 28.94515,
    "FLOW assumes air but SDS_collision_gas_mass_amu not 28.94515.")
  local mu_pa_s = FLOW.compute_mu('air', SDS_temperature_K)
  print('mu=', mu_pa_s, 'Pa s')
  FLOW.define_sphere {
    SDS = SDS,
    v_mean = 0.15,
    ri_mm = 12.7,
    ro_mm = 14.7, 
    p0_torr = SDS_pressure_torr,  -- Pressure at origin (torr)
    mu_pa_s = mu_pa_s  -- dynamic viscosity (Pa s)
  }
end
