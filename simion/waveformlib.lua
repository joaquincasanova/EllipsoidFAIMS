-- waveformlib.lua
-- Library for defining SIMION electrode potentials as segmented
-- waveforms of line segments.
--
-- Note: times are in units of microseconds.
--
-- (c) 2007-2008 Scientific Instrument Services, Inc.
-- SIMION 8.1 - PRERELEASE

-- http://simion.com/issue/490
if checkglobals then checkglobals() end

local WAVE = {}

local Type = require "simionx.Type"

-- Current waveforms to use.
-- Defined by user.
WAVE.waveforms_ = {}

-- PE surface updates are triggered every this number
-- of microseconds.  Set to nil to disable PE surface updates.
-- These are useful for display and for debugging to ensure
-- potentials are oscillating as expected.
-- Defined by user.
WAVE.pe_update_period = nil

-- Table to install segments into.
WAVE.segment = {}


-- Plots given waveforms.
-- (e.g. to Excel and/or log window).
-- This is for display/debugging purposes only.
function WAVE.plot_waveform(waves, maxtime)
  waves = waves or WAVE.waveforms_

  local EXCEL = simion.import "excellib.lua"

  local t = {}
  t.title = 'Waveform'
  t.xlabel = 'Time (usec)'
  t.ylabel = 'Potential'
  t.lines = true
  t.datasets = {}

  for iinstance,wavesi in pairs(waves) do
    for ielectrode,wave in pairs(wavesi) do
      local is_unique =
        iinstance == '*' or
        iinstance ~= '*' and
          (not waves['*'] or not waves['*'][ielectrode] or
          waves['*'][ielectrode] ~= wave)
      if is_unique then
        local dataset = {}
        dataset.header = {'Time (usec)', 'V' .. ielectrode}
        local tm = -math.huge
        while 1 do
          local lasttm = tm
          tm = wave:next_transition(tm)
          if not tm or tm == math.huge then break end
          if tm > maxtime then
            dataset[#dataset+1] = {maxtime, wave:get_potential(maxtime)}
            break
          end
          -- Allow for jumps in y by plotting left limit too.
          if tm - lasttm < math.huge then
            local tm_minus = lasttm*0.01+tm*0.99  --fudge factor-IMPROVE?
            dataset[#dataset+1] = {tm_minus, wave:get_potential(tm_minus)}
          end
          dataset[#dataset+1] = {tm, wave:get_potential(tm)}
        end
        t[#t+1] = dataset
      end
	end
  end
  
  EXCEL.plot(t)
end


-- for debugging.
local function record_potentials(t, potentials)
  local vs = ''
  for i,v in pairs(potentials) do
    vs = vs .. ',' .. v
  end
--  print('DEBUG:t=' .. t .. ',v=' .. vs)
end


-- Gets index n of piece of waveform (wave) containing time t.
--
-- That is, waves[n].time <= t < waves[n+1].time,
-- provided we conventionally define waves[0] = -infinity and
-- wave[#waves+1] = infinity.
local function get_piece(wave, t)
  local n = 0
  for m = 1, #wave do
    if wave[m].time > t then break end
    n = m
  end
  return n
end


-- SIMION segment called to override time-step size.
--
-- Ensures time-step stops precisely on end of current waveform part.
-- This is not essential (i.e. can be omitted) but does improve
-- accuracy for stepped waveforms.
local min = math.min
function WAVE.segment.tstep_adjust()
  local tm = ion_time_of_flight
  local waves = WAVE.waveforms_[ion_instance] or WAVE.waveforms_['*']
  if waves then
    for i,wave in pairs(waves) do
      local tmtransition = wave:next_transition(tm)
      if tmtransition then
        ion_time_step = min(ion_time_step, tmtransition - tm)
      end
    end
  end
end


-- SIMION segment called to adjust electrode potentials.
--
-- Updates each electrode potential.
local lastv = {}  -- store potentials for debugging only
function WAVE.segment.fast_adjust()
  local t = ion_time_of_flight
  local waves = WAVE.waveforms_[ion_instance] or WAVE.waveforms_['*']
  if waves then
    for i,wave in pairs(waves) do
      local v = wave:get_potential(t)
      if v then
        -- print('DEBUG:',t, v)
        adj_elect[i] = v   -- set potential
        lastv[i] = v
      end
    end
  end
end


-- SIMION segment called on each time step.
--
-- Updates PE surface display and print potentials.  This is done
-- periodically and also whenever transition between parts of a
-- waveform are done.  This is not essential (i.e. can be omitted) but
-- can be useful for display and debugging.
local piece_last = {}
local last_tof = -math.huge
local abs = math.abs
function WAVE.segment.other_actions()
  local period = WAVE.pe_update_period
  local is_update
  if abs(ion_time_of_flight - last_tof) > period then
    is_update = true
  else
    local waves = WAVE.waveforms_[ion_instance] or WAVE.waveforms_['*']
    if waves then
      for i,wave in pairs(waves) do
        local nt = wave:next_transition(ion_time_of_flight)
        if nt and nt >= ion_time_of_flight - ion_time_step and
                  nt < ion_time_of_flight
        then
          is_update = true
        end
      end
    end
  end
  if is_update then
    sim_update_pe_surface = 1
    last_tof = ion_time_of_flight
    if debug then
      record_potentials(ion_time_of_flight - ion_time_step, lastv)
    end
  end
end


local function merge_segments(t)
  for name,newseg in pairs(t) do
    local oldseg = segment[name]
    segment[name] =
      oldseg and function() oldseg(); newseg() end
             or  newseg
  end
end


-- Install SIMION segments.
-- This may be passed a table of name-value pairs of parameters
-- defined above: waves, pe_update_period and segment.
-- This must be called at the top-level, outside of any segments.
function WAVE.install(t)
  t = t or {}
  for k,v in pairs(t) do
    WAVE[k] = v
  end
  if WAVE[1] then
    WAVE.waveforms_ = WAVE[1]
    WAVE[1] = nil
  end

  merge_segments(WAVE.segment)
end


-- Sets waveform.  This may be called after WAVE.install to change
-- waveform.
function WAVE.set_waveform(waveforms)
  WAVE.waveforms_ = waveforms
end


-- metatable for lines objects.
local lines_mt = {}
lines_mt.__index = lines_mt

function lines_mt.__tostring() return "lines" end

-- Gets potential at time t for waveform (lines).
function lines_mt.get_potential(lines, tm)
  -- Locate current line piece [n, n+1] of the waveform.
  local n = get_piece(lines, tm)

  if n >= 1 and n < #lines then
    -- Obtain points (t1,v1) and (t2,v2) of that line segment.
    local wm,wp = lines[n],lines[n+1]
    local t1,t2, v1,v2 = wm.time,wp.time, wm.potential,wp.potential

    -- Linearly interpolate potential over the line segment.
    local v = v1 + (tm - t1) * ((v2-v1)/(t2-t1))
    return v
  elseif n == #lines and #lines > 0 and lines[n].time == tm then
    local wm = lines[n]
    return wm.potential
  end
end



-- Gets next transition time
function lines_mt.next_transition(lines, tm)
  local n = get_piece(lines, tm)
  if n < #lines then
    return lines[n+1].time  -- next point
  end
end


-- Defines a series of lines segment in a waveform.
-- Example:
--   WAVE.lines {
--     {time=0,  potential=0};
--     {time=10, potential=3};
--     {time=30, potential=-3};
--     {time=40, potential=0};
--   }
local T
function WAVE.lines(t)
  -- param check
  if not Type.is_array(t) then
    error("lines not passed array", 2)
  end
  T = T or {
    time = Type.number;
    potential = Type.number;
  }
  for i,v in ipairs(t) do
    local ok,message = Type.is_param_table(T, v)
    if not ok then
      error(message, 2)
    end
  end

  if #t > 0 and t[1].time ~= 0 then
    error("starting time of first line segment must be 0", 2)
  end

  local duration = #t > 0 and t[#t].time or 0
  t.duration = duration

  return setmetatable(t, lines_mt)
end

local argt
function WAVE.constant(t)
  argt = argt or Type {
    duration = Type.nonnegative_number + Type['nil'],
    potential = Type.number
  }
  argt:check(t)

  return WAVE.lines {
    {time=0, potential=t.potential},
    {time=t.duration or math.huge, potential=t.potential}
  }
end


-- helper function for sequence.
-- Gets index n of piece of waveform (wave) containing time t.
-- or nil if none.
-- n is the maximum value satifying the condition
--   begin_time(wave[n]) <= t <= end_time(wave[n])
local function get_piece2(wave, t)
  local wave2, toffset
  local time = 0
  for m = 1, #wave do
    if time > t then break end
    -- invar: time <= t
    local nexttime = time + wave[m].duration
    if t <= nexttime then
      wave2 = wave[m]
      toffset = time
    end
    time = nexttime
  end
  return wave2,toffset
end


-- metatable for sequence objects.
local sequence_mt = {}
sequence_mt.__index = sequence_mt

function sequence_mt.__tostring() return "sequence" end


-- Gets potential at time t for waveform (sequence).
function sequence_mt.get_potential(sequence, tm)
  -- Locate current line piece [n, n+1] of the waveform.
  local wave, toffset = get_piece2(sequence, tm)
  if wave then
    return wave:get_potential(tm - toffset)
  end
end


--FIX:check
-- Gets next transition time
function sequence_mt.next_transition(sequence, tm)
  if tm < 0 then return 0 end
  local wave, toffset = get_piece2(sequence, tm)
  if wave then
    local nt = wave:next_transition(tm - toffset)
    if nt then nt = nt + toffset end
    return nt
  end
end

-- Defines a series of sequence segment in a waveform.
-- Example:
--   WAVE.sequence {
--     .....
--   }
local T
function WAVE.sequence(t)
  -- param check
  if not Type.is_array(t) then
    error("sequence not passed array", 2)
  end

  local duration = 0
  for _,wave in ipairs(t) do
    duration = duration + wave.duration
  end
  t.duration = duration

  return setmetatable(t, sequence_mt)
end



-- metatable for loop objects.
local loop_mt = {}
loop_mt.__index = loop_mt

function loop_mt.__tostring() return "loop" end

-- Gets potential at time t.
function loop_mt.get_potential(loop, tm)
  if tm >= 0 and tm <= loop.duration then
    local wave2 = loop[1]
    local tm2 = tm % wave2.duration
    return wave2:get_potential(tm2)
  end
end

-- Gets next transition time
function loop_mt.next_transition(loop, tm)
  if tm < loop.duration then
    local n = loop.n
    local wave2 = loop[1]
    local duration2 = wave2.duration
    local shift
	tm = tm * (1+1E-14)
    local tm2
    if tm >= 0 then
      shift = math.floor(tm / duration2) * duration2
      tm2 = tm % duration2
    else
      shift = 0
      tm2 = tm
    end
    local tnext = wave2:next_transition(tm2)
    if tnext then tnext = tnext + shift end
    -- numerical roundoff error can cause tnext == tm,
    -- which results in zero time-step size (tnext - tm),
    -- and the particle freezing.  This corrects that.
    if tnext == tm then
      tnext = tnext + duration2
    end
    return tnext
  end
end


-- Defines an repetition of the contained waveforms.
-- Example:
--  WAVE.loop(10) {
--    WAVE.lines {
--      {time=0,  potential=0};
--      {time=10, potential=3};
--      {time=30, potential=-3};
--      {time=40, potential=0};
--    }
--  }
function WAVE.loop(n)
  -- param check
  if not Type.is_nonnegative_integer(n) then
    error("loop not passed non-negative integer", 2)
  end
  return function(t)
    if not Type.is_array(t) then
      error("loop not passed array", 2)
    end
    local wave2 = t[1]
    if not Type.is_table(wave2) then
      error("loop not passed array of objects", 2)
    end
    local duration = wave2.duration * n
    return setmetatable({n=n, duration=duration, wave2}, loop_mt)
  end
end

-- Defines an electrode in a waveform.
-- Example:
--    WAVE.electrode(2,3) {
--      WAVE.lines { ..... }
--    }
-- electrode(i,j) specifies electrode number i in PA instance
-- number j.  If j is omitted, it applies to all electrodes with
-- that number (in any PA instance).
function WAVE.electrode(i, j, ...)
  -- param check
  if select('#', ...) ~= 0 then
    error('extra parameters to electrode', 2)
  end
  if not Type.is_nonnegative_integer(i) then
    error('electrode number not non-negative integer', 2)
  end
  if not (Type.is_nonnegative_integer(j) or j == nil) then
    error('PA instance number not non-negative integer or nil', 2)
  end

  return function(t)
    -- param check
    if not(Type.is_array(t) or Type.is_table(t) and t.potential) then
      error('electrode not passed an array or a table with a "potential" field.', 2)
    end

    if #t == 0 and t.potential then
      t = WAVE.lines {
         {time=0,         potential=t.potential},
         {time=math.huge, potential=t.potential}
      }
    else
      t = WAVE.sequence(t)
    end

    local wave = t

    return {'electrode', i, j, wave}
  end
end

local waveforms_mt = {}

-- Defines waveforms for multiple electrodes.
-- Example:
--   WAVE.waveforms {
--     WAVE.electrode(1) { ..... };
--     WAVE.electrode(2) { ..... };
--     ...
--   }
function WAVE.waveforms(t)
  -- param check
  if not Type.is_array(t) then
    error('electrode not passed an array', 2)
  end
  for i,v in ipairs(t) do
    if type(v) ~= 'table' or v[1] ~= 'electrode' then
      error('invalid object in waveform (expected electrode)', 2)
    end
  end

  -- Build index structure used by SIMION segments.
  local result = {}
  for i,v in ipairs(t) do
    local ielectrode, iinstance, electrode = v[2], v[3], v[4]
    iinstance = iinstance or '*' -- '*' = any instance
    result[iinstance] = result[iinstance] or {}
    result[iinstance][ielectrode] = electrode
  end

  -- Propagate electrode definitions shared by all PA instances ('*')
  -- to PA instances where those electrode definitions are absent.
  for ielectrode,electrode in pairs(result['*'] or {}) do
    for iinstance,wavesi in pairs(result) do
      if not wavesi[ielectrode] then
        wavesi[ielectrode] = electrode
      end
    end
  end

  return setmetatable(result, waveforms_mt)
end


return WAVE
