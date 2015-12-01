-- excellib.lua
-- Utility library for recording and plotting data in Excel
-- (e.g. scatterplots).
--
-- Example usage (plots three points with a line between them):
--
--  local EXCEL = dofile "excellib.lua"
--  EXCEL.plot {lines=true, {1,2}, {4,5}, {8,10}}
--
-- D.Manura, 2008-05/2008-08
-- (c) 2008 Scientific Instrument Services, Inc. (Licensed under SIMION 8.0)

local EXCEL = {}


-- Excel plot types.
local xlXYScatter = -4169    -- XY scatterplot
local xlXYScatterLines = 74  -- XY scatterplot with lines


-- Returns true/false whether Excel is installed.
function EXCEL.is_excel_installed()
  return luacom.CLSIDfromProgID "Excel.Application" ~= nil
end


-- Returns an Excel object.
-- Keeps at most one Excel object open at a time (_G.excel).
-- Raises on error.
function EXCEL.get_excel()
  if not EXCEL.is_excel_installed() then
    error("Excel does not appear installed.")
  end
  -- If Excel object exists but doesn't work, perhaps the Excel process
  -- was prematurely terminated.  If so, release the Excel object.
  if _G.excel and not pcall(function() local v = _G.excel.Visible end) then
    _G.excel = nil
  end
  -- Attempt to create new Excel object if it doesn't exist.
  _G.excel = _G.excel or luacom.CreateObject("Excel.Application")
  if not _G.excel then
    error("Could not create Excel object.")
  end
  return _G.excel
end


-- Plots data from table t in Excel.
--
-- t is normally an array of rows of column values.  Alternately, t
-- may be an array of arrays of rows of column values.
-- t may optionally also contain these fields:
--
--   header - array of column headers.
--   title - title for plot.  Defaults to none if omitted.
--   xlabel - x label for plot.  Defaults to t.header[1] if omitted.
--   ylabel - y label for plot.  Defaults to t.header[2] if omitted.
--   lines - Whether to draw lines - true/false (default false)
--
-- Examples:
--
--   Plot points (1,2) and (4,5)
--   EXCEL.plot {{1,2}, {4,5}}
--
--   -- Plot with additional parameters
--   EXCEL.plot {header={'time', 'speed'},
--               title='my plot', xlabel='t', ylabel='s', lines=true,
--               {1,2}, {4,5}}
--
--   -- Plot two data series having same X values.
--   -- First column is X.  Second and third columns are Y for each series.
--   EXCEL.plot {{1,2,3}, {4,5,6}}
--
--   -- Plot two data series having possibly different X values.
--   -- Two independent sets of data:
--   -- First column is X.  Second column is Y.
--   EXCEL.plot {{{1,2},{4,5}}, {{1,3},{4,6}}}
--
function EXCEL.plot(t)
  local excel = EXCEL.get_excel()

  -- Set up worksheet.
  excel.Visible = true
  local wb = excel.Workbooks:Add()
  local ws = wb.Worksheets(1)

  -- Set up chart.
  local chart = excel.Charts:Add()
  chart.ChartType = t.lines and xlXYScatterLines or xlXYScatter

  -- Normalize table.
  local datasets = t
  if t[1] and t[1][1] and type(t[1][1]) ~= 'table' then
    datasets = {t}
  end

  -- For each dataset...
  local icol = 1
  for _,dataset in ipairs(datasets) do
    -- number of rows and columns in data.
    local nrows = #dataset
    local ncols = #(dataset.header or dataset[1])

    -- Transfer header labels and data to Excel.
    if dataset.header then
      ws:Range(ws.Cells(1,icol), ws.Cells(1, icol+ncols-1)).Value2 = dataset.header
    end
    if nrows > 0 then
      -- workaround for LuaCOM bug http://simion.com/issue/495.2
      -- is to transfer data in chunks.
      -- OLD: ws:Range(ws.Cells(2,icol), ws.Cells(nrows+1, icol+ncols-1)).Value2
      --       = dataset
      local i=1
      while i <= nrows do
        local copy = {}
        for j=1,500 do copy[j] = dataset[i + j - 1] end

--FIX:row/column order in Excel 2007?
        ws:Range(ws.Cells(1+i,1), ws.Cells(1+i+#copy-1, ncols)).Value2 = copy
        i = i + 500
      end

      -- Define chart data sources.
      for i=1,ncols-1 do
        local series = chart.SeriesCollection(chart):NewSeries()
        series.Name =    ws:Range(ws.Cells(1, icol+i), ws.Cells(1, icol+i))
        series.XValues = ws:Range(ws.Cells(2, icol),   ws.Cells(2+nrows, icol))
        series.Values  = ws:Range(ws.Cells(2, icol+i), ws.Cells(2+nrows, icol+i))
      end
    end

    icol = icol + ncols
  end

  -- Define chart options.
  if t.title then
    chart.HasTitle = true
    chart.ChartTitle:Characters().Text = t.title
  end
  local xlabel = t.xlabel or t.header and t.header[1]
  if xlabel then
    chart.Axes(1,1).HasTitle = true
    chart.Axes(1,1).AxisTitle:Characters().Text = xlabel
  end
  local ylabel = t.ylabel or t.header and t.header[2]
  if ylabel then
    chart.Axes(1,2).HasTitle = true
    chart.Axes(1,2).AxisTitle:Characters().Text = ylabel
  end

  wb.Saved = true  -- prevent asking to save.
end



return EXCEL

