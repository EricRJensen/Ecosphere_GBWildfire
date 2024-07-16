// Variables Calculated:
// -------- Temperature 
// Mean temperature: October-September (MT)
// Spring temperature: April-May (SprT)
// Summer temperature: July-August (SumT)
// Winter temperature: November-February (WinT)
// Growing season temperature: April-September (GrST)
// Mean summer maximum temperature (SumMaxT)
// Mean winter minimum temperature (WinMaxT)
// Summer-winter temperature differential (SumWinDiff)

// -------- Precipitation
// Mean annual precip: October-September (water year) (MAP)
// Spring precip: April-May (SprP)
// Summer precip: July-August (SumP)
// Winter preciptiation November-February (WinP)
// Growing season precip: April-September (GrSP)

// -------- Others:
// Frost free period: Length of the frost-free period (days) (FFree)
// Julian date of the last freezing date of spring (SprF)
// Julian date of the first freezing date of autumn (AutF)
// Mean maximum SWE (SWE)
// Degree-days >5 degrees C (DD5)
// Degree-days <0 degrees C (DD0)
// Degree-days >5 degrees C accumulating within the frost-free period (DD5ff)
// Degree-days <0 degress C min temp (DD0Min)
// Annual Dryness index: dd5/map or sqrt(dd5)/map (once named ami annual moisture index) (ADI)
// Summer Dryness index: gsdd5/gsp or sqrt(gsdd5)/gsp (once named smi, summer moisture index) (SDI)
// Mean maximum SWE (SWE)

// Link to the Rehfeldt model document: https://www.fs.fed.us/rm/pubs/rmrs_gtr165.pdf



// -------------------------  Import and Preprocess layers ----------------------------

// Import datasets
var Daymet = ee.ImageCollection("NASA/ORNL/DAYMET_V3"),
    GBbounds = ee.FeatureCollection("users/ericjensen41_default/Thesis/Project_Boundaries/GBbounds"),
    AIMplots = ee.FeatureCollection("users/ericjensen41_default/Thesis/Plots/Allplots"),
    fire_pts = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/Fire_points"),
    SaddleDraw = ee.FeatureCollection("users/ericjensen41_default/Thesis/Project_Boundaries/SaddleDraw");

// Select fire to export climate averages for 
var SaddleDraw = SaddleDraw.first()
var GBbb = SaddleDraw.geometry()
print(GBbb)


// ----------------- Preprocess collection for making calcs ----------------------

// Add a mean temperature band
var calc_tmean = function(img){
      var tmean = img.expression(
        '(tmax + tmin) / 2', {
          'tmax': img.select('tmax'),
          'tmin': img.select('tmin')}).rename('tmean')
      var imgSmean = img.addBands(tmean)
      return(imgSmean)
}
var Daymet = Daymet.map(calc_tmean)

// Add a binary band for frost free pixels for each date
//  1 = frost-free, 0 = fros
var calc_frostfree = function(img){
  var frostfree = img.select('tmin').gte(0).rename('frostfree')
  var img_wFF = img.addBands(frostfree)
  return(img_wFF)
}
var Daymet = Daymet.map(calc_frostfree)

// Add day of year band
var addDate = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().rename('doy')
  doyBand = doyBand.updateMask(image.select('prcp').mask())
  
  return image.addBands(doyBand);
};
var Daymet = Daymet.map(addDate)


// -------------------------- Temporal sequence ---------------------------------

// Filter to 30 year period
var startDate = '1991-01-01'
var endDate = '2020-12-31'
var DM = Daymet.filterDate(startDate,endDate).filterBounds(GBbb).select('prcp', 'tmean', 'swe', 'tmin', 'tmax', 'frostfree', 'doy')

// Define start and end years for calculation
var startYear = 1990;
var endYear = 2019;

// Make a list of years to generate composites for.
var yearList = ee.List.sequence(startYear, endYear);


// ----------------------- Seasonal calendar filters ---------------------------

// Months match Rehfeldt 2006, link above

var filterSpr = ee.Filter.calendarRange(4, 5, 'month') // Spring filter
var filterSum = ee.Filter.calendarRange(7, 8, 'month') // Summer filter
var filterWin = ee.Filter.calendarRange(11, 2, 'month') // Winter filter
var filterGrS = ee.Filter.calendarRange(4, 9, 'month') // Growing filter
var filterSpF = ee.Filter.calendarRange(1, 7, 'month') // Last spring frost filter // January - July 
var filterAuF = ee.Filter.calendarRange(8, 12, 'month') // First autumn frost filter // August - December


// ----------------------- Calculate 30-year means ----------------------------

// ---------------------- Calculate precipitation -----------------------------
// Annual precipitation (MAP)
var MAP = DM.select('prcp').reduce(ee.Reducer.sum()).divide(30).rename('MAP').clip(GBbb)

// Spring precipitation (SprP)
var Spr_IC = DM.filter(filterSpr)
var SprP = Spr_IC.select('prcp').reduce(ee.Reducer.sum()).divide(30).rename('SprP').clip(GBbb)

// Summer precipitation (SumP)
var Sum_IC = DM.filter(filterSum)
var SumP = Sum_IC.select('prcp').reduce(ee.Reducer.sum()).divide(30).rename('SumP').clip(GBbb)

// Growing season precipitation (GrSP)
var GrS_IC = DM.filter(filterGrS)
var GrSP = GrS_IC.select('prcp').reduce(ee.Reducer.sum()).divide(30).rename('GrSP').clip(GBbb)

// Winter precipitation (WinP)
var Win_IC = DM.filter(filterWin)
var WinP = Win_IC.select('prcp').reduce(ee.Reducer.sum()).divide(30).rename('WinP').clip(GBbb)


// ------------ Calculate temperature ---------------------
//  Mean temperature (MT)
var MT = DM.select('tmean').reduce(ee.Reducer.mean()).rename('MT').clip(GBbb)

// Mean spring temperature (SprT)
var SprT = Spr_IC.select('tmean').reduce(ee.Reducer.mean()).rename('SprT').clip(GBbb)

// Mean summer temperature (SumT)
var SumT = Sum_IC.select('tmean').reduce(ee.Reducer.mean()).rename('SumT').clip(GBbb)

// Growing season temperature (GrST)
var GrST = GrS_IC.select('tmean').reduce(ee.Reducer.mean()).rename('GrST').clip(GBbb)

// Mean winter temperature (WinT)
var WinT = Win_IC.select('tmean').reduce(ee.Reducer.mean()).rename('WinT').clip(GBbb)

// Summer-winter temperature differential (SumWinDiff)
var SumWinDiff = ee.Image(SumT.subtract(WinT)).rename('SumWinDiff')


// ------------ Calculate seasonal extreme temperatures -------------------
// ##### Mean summer maximum (SumMaxT) #####
// Reducer for image collection to annual maximum
var SumMaxTReducer = ee.Reducer.max();

// Map reducer over the list of years to generate image collections of annual maximum images from November (previous year) to April (current year)
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var SumMaxTCol = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-06-01' // Start in November
  var enddate = '-10-01'// End in April
  var yearCol = DM.filterDate(ee.String(year_num.toInt()).cat(startdate), ee.String(year_num.toInt()).cat(enddate)).select("tmax");
  var yearComp = yearCol.reduce(SumMaxTReducer).clip(GBbb);
  return yearComp.set({
    'year': year});
}));

// Create image of mean summer maximum temperatures
var SumMaxT = SumMaxTCol.reduce(ee.Reducer.mean()).rename('SumMaxT')
// print(SumMaxT)
// Map.addLayer(SumMaxT, {min:25, max:45}, 'SumMax')


// ##### Mean winter minimum (WinMaxT) #####
// Reducer for image collection to annual minimum
var WinMinTReducer = ee.Reducer.min();

// Map reducer over the list of years to generate image collections of annual maximum images from November (previous year) to April (current year)
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var WinMinTCol = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-11-01' // Start in November
  var enddate = '-04-01'// End in April
  var yearCol = DM.filterDate(ee.String(year_num.subtract(1).toInt()).cat(startdate), ee.String(year_num.toInt()).cat(enddate)).select("tmin");
  var yearComp = yearCol.reduce(WinMinTReducer).clip(GBbb);
  return yearComp.set({
    'year': year});
}));

// Create image of mean winter minimum temperatures
var WinMinT = WinMinTCol.reduce(ee.Reducer.mean()).rename('WinMinT')


// ------------ Calculate frost free period, last frost and first frost -------------------
// SF: rehfeldt 2006 has this crazy way of estimating ffp using stepwise regression because
//     they don't have daily maps. Since we have daily maps this way looks good.
//     In reading through the code I didn't see any problems but I didn't do any thorough checks or anything.

// Mask pixels with 'tmin' < 0 and keep only pixels when frost occurred and only the DOY band
var DM_frost  = DM.map(function(image){
  return image.updateMask(image.select('tmin').lt(0)).select('doy');
});

// Filter by date to create spring frost and autumn frost collections
var SprFCol = DM_frost.filter(filterSprF)
var AutFCol = DM_frost.filter(filterAutF)

// Define reducers for last spring frost and first autumn frost
var SprFReducer = ee.Reducer.lastNonNull();
var AutFReducer = ee.Reducer.firstNonNull();

// Map reducer over the list of years to generate image collections for spring and autumn frost with a composite for each year.
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var SprFCol = ee.ImageCollection(yearList.map(function(year){
  var yearCol = SprFCol.filter(ee.Filter.calendarRange(year, year, 'year'));
  var yearComp = yearCol.reduce(SprFReducer).clip(GBbb);
  return yearComp.set({
    'year': year});
}));

var AutFCol = ee.ImageCollection(yearList.map(function(year){
  var yearCol = AutFCol.filter(ee.Filter.calendarRange(year, year, 'year'));
  var yearComp = yearCol.reduce(AutFReducer).clip(GBbb);
  return yearComp.set({
    'year': year});
}));

// Reduce Spring and Autumn frost image collections to produce the 30-year mean image
var SprF = SprFCol.reduce(ee.Reducer.mean()).rename('SprF')
var AutF = AutFCol.reduce(ee.Reducer.mean()).rename('AutF')
// Calculate mean frost free period to produce the 30-year mean image
var FFree = AutF.subtract(SprF).rename('FFree')


// ------------- Calculate growing degree days ------------------------

// ------------- Greater than 5-degree days within the frost free period
// Function to filter and calculate 5-degree-days for each day during frost free period
var calc_dd5_ff = function(img, year ){
  //  Add degree day band
  var dd5 = ee.Image(img.select('tmean').subtract(5)).rename(['dd5ff']) // 5-degree-day calculation
  var img_dd5 = img.addBands(dd5).select(['doy', 'dd5ff']) // add 5-degree-day band to image 
  
  // Filter dd5 values lt 0
  var dd5_gte0_mask = img_dd5.select('dd5ff').gte(0) // Mask pixels of dd5 less than 0
  var img_dd5_clean = img_dd5.updateMask(dd5_gte0_mask) // Apply mask
  
  // Using same year as the image to get ffp—SF
  var SprF_mask = img_dd5_clean.select('doy').gt(SprFCol.filterMetadata('year', 'equals', year).first()) // filter pixels based on doy < doy of last spring frost
  var AutF_mask = img_dd5_clean.select('doy').lt(AutFCol.filterMetadata('year', 'equals', year).first()) 
  
  var img_dd5_ff = img_dd5_clean.updateMask(SprF_mask).updateMask(AutF_mask) // apply masks to image to create frostfree image
  
  return img_dd5_ff
}

// Reducer for image collection to annual sum
var DD5Reducer = ee.Reducer.sum();

// Map degree day filtering function and sum reducer over annual image collection to produce annual accumulations of 5-degree-days during frost free period
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var DD5ffCol = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-01-01' // Entire year
  var enddate = '-01-01'// Entire year
  var yearCol = DM.filterDate(ee.String(year_num.toInt()).cat(startdate), ee.String(year_num.add(1).toInt()).cat(enddate)); // filter by dates
  var DD5Col = yearCol.map(function(i){return calc_dd5_ff(i, year)}) // Map DD5 calculation function over annual image collection // SF: Changed to use FFP for same year
  var yearComp = DD5Col.reduce(DD5Reducer).clip(GBbb) // summing reducer to return image of accumulated degree days
  return yearComp.set({
    'year': year});
}));

// Calculate 30-year mean
var DD5ff = DD5ffCol.reduce(ee.Reducer.mean()).select('dd5ff_sum_mean').rename('DD5FF')

// print(DD5ff)
// Map.addLayer(DD5ff, {min: 300, max: 2000}, 'dd5ff')


// ------------- Greater than 5-degree days, no constraint on frost free period
// Function to filter and calculate 5-degree-days for each day
var calc_dd5 = function(img){
  //  Add degree day band
  var dd5 = ee.Image(img.select('tmean').subtract(5)).rename(['dd5']) // 5-degree-day calculation
  var img_dd5 = img.addBands(dd5).select(['doy', 'dd5']) // add 5-degree-day band to image 
  
  // Filter dd5 values lt 0
  var dd5_gte0_mask = img_dd5.select('dd5').gte(0) // Mask pixels of dd5 less than 0
  var img_dd5_clean = img_dd5.updateMask(dd5_gte0_mask) // Apply mask
  
  return img_dd5_clean
}

// Reducer for image collection to annual sum
var DD5Reducer = ee.Reducer.sum();

// Map degree day filtering function and sum reducer over annual image collection to produce annual accumulations of 5-degree-days
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var DD5Col = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-01-01' // Entire year
  var enddate = '-01-01'// Entire year
  var yearCol = DM.filterDate(ee.String(year_num.toInt()).cat(startdate), ee.String(year_num.add(1).toInt()).cat(enddate)); // filter by dates
  var DD5Col = yearCol.map(calc_dd5) // Map DD5 calculation function over annual image collection
  var yearComp = DD5Col.reduce(DD5Reducer).clip(GBbb) // summing reducer to return image of accumulated degree days
  return yearComp.set({
    'year': year});
}));

// Calculate 30-year mean
var DD5 = DD5Col.reduce(ee.Reducer.mean()).select('dd5_sum_mean').rename('DD5')

// print(DD5)
// Map.addLayer(DD5, {min: 300, max: 2000}, 'dd5')


// ------------- Less than 0-degree days, no constraint on frost-free period
// Function to filter and calculate 0-degree-days for each day
var calc_dd0 = function(img){
  //  Add degree day band
  var dd0 = ee.Image(img.select('tmean').subtract(0)).rename(['dd0']) // 0-degree-day calculation - no calculation needed
  var img_dd0 = img.addBands(dd0).select(['doy', 'dd0']) // add 0-degree-day band to image 
  
  // Filter dd0 values gt 0
  var dd0_lte0_mask = img_dd0.select('dd0').lte(0) // Mask pixels of dd0 greater than 0
  var img_dd0_clean = img_dd0.updateMask(dd0_lte0_mask) // Apply mask
  
  // Convert negative values to absolute values (degree days)
  // SF: Why absolute? 
  // EJ: In the calculation of degree days in Rehfeldt it says to subtract the threshold temperature (0 degrees) from the mean temperature 
  // for that day. So, if the mean temperature is -5 °C then the calculation is -5 - 0 = -5. However, my understanding is that for degree-days 
  // below 0°C (seemingly synonymously known as freezing degree days) those values should be expressed as a positive value So for that day 
  // described previously there was an accumulation of 5 freezing degree days. If you see the graph on page 17 in Rehfeldt degree-days <0°C are 
  // expressed positively. This seems to be backed up elsewhere: https://sites.google.com/site/cryospherecomputing/fdd
  var img_dd0_abs = img_dd0_clean.abs()
  
  return img_dd0_abs
}

// Reducer for image collection to annual sum
var DD0Reducer = ee.Reducer.sum();

// Map degree day filtering function and sum reducer over annual image collection to produce annual accumulations of 5-degree-days
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var DD0Col = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-08-01' // Entire year from August to August: https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-14-0119.1
  var enddate = '-08-01'// See above
  var yearCol = DM.filterDate(ee.String(year_num.subtract(1).toInt()).cat(startdate), ee.String(year_num.toInt()).cat(enddate)); // filter by dates
  var DD0Col = yearCol.map(calc_dd0) // Map DD5 calculation function over annual image collection
  var yearComp = DD0Col.reduce(DD0Reducer).clip(GBbb) // summing reducer to return image of accumulated degree days
  return yearComp.set({
    'year': year});
}));

// Calculate 30-year mean
var DD0 = DD0Col.reduce(ee.Reducer.mean()).select('dd0_sum_mean').rename('DD0')

// print(DD0)
// Map.addLayer(DD0, {min: 0, max: 800}, 'dd0')


// ------------- Less than 0-degree days for minimum temperature, no constraint on frost-free period
// Function to filter and calculate 0-degree-days for each day
var calc_dd0_min = function(img){
  //  Add degree day band
  var dd0 = ee.Image(img.select('tmin').subtract(0)).rename(['dd0']) // 0-degree-day calculation - no calculation needed
  var img_dd0 = img.addBands(dd0).select(['doy', 'dd0']) // add 0-degree-day band to image 
  
  // Filter dd0 values gt 0
  var dd0_lte0_mask = img_dd0.select('dd0').lte(0) // Mask pixels of dd0 greater than 0
  var img_dd0_clean = img_dd0.updateMask(dd0_lte0_mask) // Apply mask
  
  // Convert negative values to absolute values (degree days)
  // SF: Why absolute? 
  // EJ: In the calculation of degree days in Rehfeldt it says to subtract the threshold temperature (0 degrees) from the mean temperature 
  // for that day. So, if the mean temperature is -5 °C then the calculation is -5 - 0 = -5. However, my understanding is that for degree-days 
  // below 0°C (seemingly synonymously known as freezing degree days) those values should be expressed as a positive value So for that day 
  // described previously there was an accumulation of 5 freezing degree days. If you see the graph on page 17 in Rehfeldt degree-days <0°C are 
  // expressed positively. This seems to be backed up elsewhere: https://sites.google.com/site/cryospherecomputing/fdd
  var img_dd0_abs = img_dd0_clean.abs()
  
  return img_dd0_abs
}

// Reducer for image collection to annual sum
var DD0Reducer = ee.Reducer.sum();

// Map degree day filtering function and sum reducer over annual image collection to produce annual accumulations of 5-degree-days
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var DD0MinCol = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-08-01' // Entire year from August to August: https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-14-0119.1
  var enddate = '-08-01'// See above
  var yearCol = DM.filterDate(ee.String(year_num.subtract(1).toInt()).cat(startdate), ee.String(year_num.toInt()).cat(enddate)); // filter by dates
  var DD0Col = yearCol.map(calc_dd0_min) // Map DD5 calculation function over annual image collection
  var yearComp = DD0Col.reduce(DD0Reducer).clip(GBbb) // summing reducer to return image of accumulated degree days
  return yearComp.set({
    'year': year});
}));

// Calculate 30-year mean
var DD0Min = DD0MinCol.reduce(ee.Reducer.mean()).select('dd0_sum_mean').rename('DD0Min')

// print(DD0)
// Map.addLayer(DD0, {min: 0, max: 800}, 'dd0')


// ------------------------- Calculate dryness indices -----------------------------
// ##### Annual dryness index #####
var adi = ee.Image(DD5.divide(MAP)).rename('ADI')

// ##### Summer dryness index #####
var sdi = ee.Image(DD5ff.divide(GrSP)).rename('SDI')

// Map.addLayer(adi, {min:0, max:30}, 'adi')
// Map.addLayer(sdi, {min:0, max:30}, 'sdi')


// ---------------- Calculate mean maximum snow water equivalent -------------------
// Reducer for image collection to annual maximum
var SWEReducer = ee.Reducer.max();

// Map reducer over the list of years to generate image collections of annual maximum images from November (previous year) to April (current year)
// Stack exchange: https://gis.stackexchange.com/questions/340433/making-intra-annual-image-composites-for-a-series-of-years-in-google-earth-engin
var SWECol = ee.ImageCollection(yearList.map(function(year){
  var year_num = ee.Number(year)
  var startdate = '-11-01' // Start in November
  var enddate = '-04-01'// End in April
  var yearCol = DM.filterDate(ee.String(year_num.subtract(1).toInt()).cat(startdate), ee.String(year_num.toInt()).cat(startdate)).select("swe");
  var yearComp = yearCol.reduce(SWEReducer).clip(GBbb);
  return yearComp.set({
    'year': year});
}));

var SWE = SWECol.reduce(ee.Reducer.mean()).rename('SWEMax')


// ------------ Compile climate means images into single image for reducing regions and exporting over -------------

// List of climate means images
var CliMeans_list = ee.List([MAP, SprP, SumP, GrSP, WinP, MT, SprT, SumT, GrST, WinT, SumWinDiff, SumMaxT, WinMinT, SprF, AutF, FFree, 
DD5ff, DD5, DD0, DD0Min, adi, sdi, SWE])

var CliNames_list = ['MAP', 'SprP', 'SumP', 'GrSP', 'WinP', 'MT', 'SprT', 'SumT', 'GrST', 'WinT', 'SumWinDiff', 'SumMaxT', 'WinMinT', 'SprF', 'AutF', 'FFree', 
'DD5ff', 'DD5', 'DD0', 'DD0Min', 'ADI', 'SDI', 'SWE']

// print(CliMeans_IC)
// print(CliMeans_IC.first().bandNames().get(0))


// Client-side for loop to produce exports; only looping over client-side JS arrays
for (var i = 0; i < 23; i++) {
  // Get image
  var img = ee.Image(CliMeans_list.get(i))
  
  // Run reduce regions on image
  var AIM_rr = img.reduceRegions({
    collection: fire_pts,
    reducer: ee.Reducer.mean(),
    scale: 30,
    tileScale: 16})
  
  // Export results to drive
  Export.table.toDrive({
      collection: AIM_rr,
      description: CliNames_list[i] + '_fires',
      fileFormat: 'csv'}) 
}


// -------------------------------- EXPORT FOR R&R MODELING --------------------------------

// Create single image to export by
var CliMeans_I = ee.ImageCollection([adi, DD0, DD5, SprP, SumP, WinP, WinT])
                        .toBands()

// // Reproject and resample image by nearest neighbors to match Landsat
var LS_ref = ee.Image('users/zackrwerner/landsat_harm_reference')
var LS_proj = LS_ref.projection() // get projection of landsat harmonized image
var CliMeans_I_resample = CliMeans_I.reproject(LS_proj) //reproject defaults to nn 

print(CliMeans_I_resample)

// // Calculate focal means at 100 meter scale using reduceNeighborhood
// var CliMeans_I_resample_focal = CliMeans_I_resample.reduceNeighborhood({
//   reducer: ee.Reducer.mean(),
//   kernel: ee.Kernel.circle(100, 'meters')})
//   .rename(['mn_adi', 'mn_SumP', 'mn_SprP', 'mn_SWE', 'mn_DD0Min', 'mn_MAP'])

// Get list of bandnames to export over
var namelist = CliMeans_I_resample.bandNames().getInfo();
print(namelist,'namelist');

for (var i = 0; i < 7; i++) {
  // Get image
  var img = CliMeans_I_resample.select(namelist[i]).clip(GBbb);
  
  // Export results to drive
  Export.image.toDrive({
      image: img,
      description: namelist[i].substr(2,10) + '_avg',
      scale: 30,
      maxPixels: 1e13,
      region: GBbb});}  


// --------------------------------- Visualize data ------------------------------------
//  Visualization parameters for temperature data
var pParams = {
  min: 100,
  max: 300,
  palette: ['8B4513','ffffff','00ff00']
};
var tParams = {
  min: 0,
  max: 20,
  palette: ['113EDA','ffffff','CF523F']
};

var ffParams = {
  min: -300,
  max: -175,
  palette: ['113EDA','ffffff','CF523F']
};
// Map.addLayer(GBbounds, {}, 'GBbounds',false)
// Map.addLayer(MAP, pParams, 'MAP',false)
// Map.addLayer(MT, tParams, 'MAM',false)
// Map.addLayer(AIMplots, {}, 'AIM',false)
// Map.addLayer(diffDays, ffParams, 'ff',false)