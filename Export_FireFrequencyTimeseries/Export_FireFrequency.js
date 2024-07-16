

// --------------------- Read-in datasets ------------------------

// Read in datasets
var AIMplots = ee.FeatureCollection("users/ericjensen41_default/Thesis/Plots/Allplots")
var Aerial = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/TreatmentPolygons/LTDL_AerialSeed")
var Drill = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/TreatmentPolygons/LTDL_DrillSeed")
var Fires = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/FirePolygons/Fires_BurnedOnce") // Select Fires frequency strata and export accordingly
// var Fires = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/FirePolygons/Fires_Unburned") // Select Fires frequency strata and export accordingly
// var Fires = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/FirePolygons/Fires_BurnedMTOnce") // Select Fires frequency strata and export accordingly
var gb = ee.FeatureCollection("users/ericjensen41_default/Thesis/Project_Boundaries/GBbounds")
var nlcd = ee.ImageCollection("USGS/NLCD")

// Read in geometries
var jensen = AIMplots.filterMetadata('Program', 'equals', 'Jensen') 
var Seeded = Drill.merge(Aerial)


// ---- Create NLCD if-ever mask of agriculture and developed ----

var gb_nlcd = nlcd.filterBounds(gb).select('landcover').toBands().clip(gb)

var const1 = ee.Image(1)
var dev_mask21 = gb_nlcd.eq(21).reduce(ee.Reducer.max())
var mask21 = const1.updateMask(dev_mask21)
var dev_mask22 = gb_nlcd.eq(22).reduce(ee.Reducer.max())
var mask22 = const1.updateMask(dev_mask22)
var dev_mask23 = gb_nlcd.eq(23).reduce(ee.Reducer.max())
var mask23 = const1.updateMask(dev_mask23)
var dev_mask24 = gb_nlcd.eq(24).reduce(ee.Reducer.max())
var mask24 = const1.updateMask(dev_mask24)
var agr_mask81 = gb_nlcd.eq(81).reduce(ee.Reducer.max())
var mask81 = const1.updateMask(agr_mask81)
var agr_mask82 = gb_nlcd.eq(82).reduce(ee.Reducer.max())
var mask82 = const1.updateMask(agr_mask82)

var nlcd_inv_mask = ee.ImageCollection([mask21,mask22,mask23,mask24,mask81,mask82]).mosaic()
var toMask = const1.updateMask(nlcd_inv_mask).unmask().remap([0,1],[1,0]).gte(1)
var nlcd_mask = const1.updateMask(toMask)


//  ---------------- Preprocess RAP Indicators -----------------

var RAP = ee.ImageCollection("projects/rangeland-analysis-platform/vegetation-cover-v2")
  // .select(['AFGC', 'PFGC', 'SHR', 'TREE']);

// Calculate normalized differenced index and add band
// Perennial - annual / perennial + annual // Note: Excluding trees
var calc_normdif = function(img){
  var per = img.select('PFGC').add(img.select('SHR')).rename('PER') // Add band for perennials
  var img_wPer = img.addBands(per)
  var normdif = img_wPer.normalizedDifference(['PER', 'AFGC']).multiply(100).rename('NDPDI') // calculate normalized difference perennial dominance index
  var img_yr = img.get('year')
  var yearband = ee.Image.constant(img_yr).toFloat().rename("Year")
  return(img_wPer.addBands(normdif).addBands(yearband))
}
var RAP = RAP.map(calc_normdif).select(['AFGC', 'BG', 'PFGC', 'SHR', 'NDPDI', 'PER', 'Year'])


// -------------------------------------------------------------
// -------------- Investigate Specific fires -------------------

// ------------ Generate Pre- and Post-fire layers -------------

// Fire names are:
// 'TUANA COMPLEX' (1995)
// 'BILK CREEK COMPLEX (DOUBLE H)' (2000)
// 'SHIRTTAIL' (1999)
// 'SADLER COMPLEX' (1999)
// 'BIG JUNIPER' (2001)
// 'SOMBRERO' (1999)
// 'WEST BASIN' (2000)

// Select the desired fire and extact objects for script to run
var Fire = Fires.filter(ee.Filter.eq('Fire_Name', 'BILK CREEK COMPLEX (DOUBLE H)'))
print(Fire)
var FireGeom = Fire.first().geometry()
var FireYear = ee.Number(Fire.first().get('Fire_Year'))
var FireYearStr = ee.String(Fire.first().get('Fire_Year'))
print('The year of the fire was', FireYear)

// Calculate pre-fire and post-fire images and recovery images
var NDPDI_pre = RAP.select('NDPDI')
            .filterDate(ee.String(FireYear.subtract(4)).cat('-01-01'), ee.String(FireYear.subtract(1)).cat('-01-01')) // Three years prior to fire // no pixels burned in three years leading up to any fires
            .toBands()
            .reduce(ee.Reducer.median()) 
var NDPDI_post = RAP.select('NDPDI').filterDate('2015-01-01', '2020-01-01').toBands().reduce(ee.Reducer.median()) // Five most recent years
var NDPDI_recov = NDPDI_post.divide(NDPDI_pre) // Percent recovery of NDPDI

// Generate recovery slopes for the fire and summarize slopes by treatment polygons
var postfire = RAP.filterMetadata('year', 'greater_than', ee.Number(FireYear))
var NDPDI_slope = postfire.select(['Year', 'NDPDI']).reduce(ee.Reducer.linearFit()).select('scale')


// -------- Calculate time-series plots of recovery for entire RAP collection -----------

// Prep unseeded image collection
var maskInside = function(image) {
  var mask = ee.Image.constant(1).clip(Seeded).mask().not()
  return image.updateMask(mask)
}
var unseeded_IC = RAP.map(maskInside)

// Clip to Fire Geometry
var clipToFire = function(img){
return(img.clip(FireGeom)) 
}

var drill_IC = RAP.map(clipToFire)
var aerial_IC = RAP.map(clipToFire)
var unseeded_IC = unseeded_IC.map(clipToFire)

// Clip drill and aerial to respective feature collections
var aerial_IC = aerial_IC.map(function(img){
  return img.clip(Aerial)
})
var drill_IC = drill_IC.map(function(img){
  return img.clip(Drill)
})

print(ui.Chart.image.series(drill_IC.select(['NDPDI', 'AFGC', 'PER']), FireGeom, ee.Reducer.mean(), 30), 'Drill seeded')
print(ui.Chart.image.series(aerial_IC.select(['NDPDI', 'AFGC', 'PER']), FireGeom, ee.Reducer.mean(), 30), 'Aerial seeded')
print(ui.Chart.image.series(unseeded_IC.select(['NDPDI', 'AFGC', 'PER']), FireGeom, ee.Reducer.mean(), 30), 'Unseeded')


// ------------- Generate histograms of NDPDI slopes along seeding types --------------------------
// Get the seeding types
var drill = NDPDI_slope.clip(Drill)
var aerial = NDPDI_slope.clip(Aerial)
var unseeded = maskInside(NDPDI_slope, Seeded).clip(FireGeom)


// Print the charts
// Visualization parameters for charts
  var options = {
    title: 'NDPDI Slopes',
    hAxis: {
      viewWindowMode: 'explicit',
      viewWindow: {
        min:-5,
        max: 5    }}};

print(ui.Chart.image.histogram({image:drill, region:FireGeom, scale:30}).setOptions(options), "NDPDI_Drill");
print(ui.Chart.image.histogram({image:aerial, region:FireGeom, scale:30}).setOptions(options), "NDPDI_Aerial");
print(ui.Chart.image.histogram({image:unseeded, region:FireGeom, scale:30}).setOptions(options), "NDPDI_Unseeded");


// ---------- Visualize the Bilk Fire analysis ---------------
// Map different recovery characteristics
Map.addLayer(NDPDI_pre.clip(FireGeom), {min:-100, max:100, palette: ['brown', 'white', 'green']}, 'NDPDI Prefire')
Map.addLayer(NDPDI_post.clip(FireGeom), {min:-100, max:100, palette: ['brown', 'white', 'green']}, 'NDPDI Postfire',false)
Map.addLayer(NDPDI_recov.clip(FireGeom), {min:-5, max:5, palette: ['brown', 'white', 'green']}, 'NDPDI recov',false)
Map.addLayer(NDPDI_slope.clip(FireGeom), {min:-5, max:5, palette: ['brown', 'white', 'green']}, 'NDPDI slope',false)

// Treatment polygons
Map.addLayer(Aerial, {color: 'green'}, 'Aerial seeding',false)
Map.addLayer(Drill, {color: 'blue'}, 'Drill seeding',false)

// View AIM plots in the Bilk Creek Fire
Map.addLayer(AIMplots.filterBounds(Fires), {}, 'Fire AIM plots', false)

// Fire data
Map.addLayer(Fire, {}, 'Fire Polygon')


// ------------------- Export CSVs of RAP indicators --------------------

// Export CSVs of NDPDI, perennials, and annuals for Jensen plots
var NDPDI_rr = RAP.select('NDPDI').toBands().reduceRegions({
  collection: jensen,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

var AFGC_rr = RAP.select('AFGC').toBands().reduceRegions({
  collection: jensen,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

var PER_rr = RAP.select('PER').toBands().reduceRegions({
  collection: jensen,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

Export.table.toDrive({
  collection: NDPDI_rr,
  description:'NDPDI_rr',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: PER_rr,
  description:'PER_rr',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: AFGC_rr,
  description:'AFGC_rr',
  fileFormat: 'CSV'
});


//  ------------------- Export CSVs for each fire ---------------------

// RAP object for exporting all pixels
var RAP_forRR = RAP.select(['AFGC', 'BG', 'PFGC', 'SHR', 'NDPDI', 'PER']).toBands().updateMask(nlcd_mask)

var fires_rr = RAP_forRR.reduceRegions({
  collection: Fires,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

print(fires_rr)

Export.table.toDrive({
  collection: fires_rr,
  description:'MTOnce_RR',
  fileFormat: 'CSV'
});

// ----------------- Export data for all fires combined --------------

// Create flat feature collection of fires for export
var Fires_flat = ee.FeatureCollection(Fires.geometry())
var MTOnce_all_rr = RAP_forRR.reduceRegions({
  collection: Fires_flat,
  reducer: ee.Reducer.mean(),
  scale: 30
});

print(MTOnce_all_rr, 'MTOnce_all_rr')

Export.table.toDrive({
  collection: MTOnce_all_rr,
  description:'MTOnce_all_rr',
  fileFormat: 'CSV'
});


// ---------------- Export seeding type for fires ----------------------

// Prep unseeded image collection
var maskInside = function(image) {
  var mask = ee.Image.constant(1).clip(Seeded).mask().not()
  return image.updateMask(mask)
}
var unseeded_IC = RAP.map(maskInside)

// Clip drill and aerial to respective feature collections
var aerial_IC = RAP.map(function(img){
  return img.clip(Aerial)
})
var drill_IC = RAP.map(function(img){
  return img.clip(Drill)
})

// Export aerial, drill, and unseeded CSVs
var drill_I = drill_IC.select(['AFGC', 'BG', 'PFGC', 'SHR', 'NDPDI', 'PER']).toBands()
var aerial_I = aerial_IC.select(['AFGC', 'BG', 'PFGC', 'SHR', 'NDPDI', 'PER']).toBands()
var unseed_I = unseeded_IC.select(['AFGC', 'BG', 'PFGC', 'SHR', 'NDPDI', 'PER']).toBands()

var drill_rr = drill_I.reduceRegions({
  collection: Fires,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

var aerial_rr = aerial_I.reduceRegions({
  collection: Fires,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

var unseed_rr = unseed_I.reduceRegions({
  collection: Fires,
  reducer: ee.Reducer.mean(),
  scale: 30,
});

Export.table.toDrive({
  collection: drill_rr,
  description:'Drill_RR',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: aerial_rr,
  description:'Aerial_RR',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: unseed_rr,
  description:'Unseeded_RR',
  fileFormat: 'CSV'
});