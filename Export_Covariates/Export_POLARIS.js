// --------------------------- Import datasets ------------------------------

// Import POLARIS collections
var BD_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_BD_00_05"),
    BD_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_BD_05_15"),
    BD_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_BD_15_30"),
    BD_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_BD_30_60"),
    Clay_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Clay_00_05"),
    Clay_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Clay_05_15"),
    Clay_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Clay_15_30"),
    Clay_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Clay_30_60"),
    OM_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_OM_00_05"),
    OM_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_OM_05_15"),
    OM_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_OM_15_30"),
    OM_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_OM_30_60"),
    Sand_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Sand_00_05"),
    Sand_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Sand_05_15"),
    Sand_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Sand_15_30"),
    Sand_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Sand_30_60"),
    Silt_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Silt_00_05"),
    Silt_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Silt_05_15"),
    Silt_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Silt_15_30"),
    Silt_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_Silt_30_60"),
    ThetaR_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaR_00_05"),
    ThetaR_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaR_05_15"),
    ThetaR_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaR_15_30"),
    ThetaR_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaR_30_60"),
    ThetaS_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaS_00_05"),
    ThetaS_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaS_05_15"),
    ThetaS_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaS_15_30"),
    ThetaS_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_ThetaS_30_60"),
    PH_00_05 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_PH_00_05"),
    PH_05_15 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_PH_05_15"),
    PH_15_30 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_PH_15_30"),
    PH_30_60 = ee.Image("users/ericjensen41_default/Thesis/Covariates/Covar_POLARIS/POLARIS_PH_30_60")

// Import other datasets
var GBbounds = ee.FeatureCollection("users/ericjensen41_default/Thesis/Project_Boundaries/GBbounds")
var fire_pts = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/Fire_points")
var fire_pts_all = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/Fire_pointsAll")

// // Buffer function to apply to fire points
// var buffer100 = function(f){
//     return f.buffer(100) }
// fire_polys = fire_pts.map(buffer100)

// Create image collection from POLARIS images
var POLARIS_IC = ee.ImageCollection([BD_00_05, BD_05_15, BD_15_30, BD_30_60, Clay_00_05, Clay_05_15, Clay_15_30, Clay_30_60,
                  OM_00_05, OM_05_15, OM_15_30, OM_30_60, Sand_00_05, Sand_05_15, Sand_15_30, Sand_30_60, Silt_00_05, Silt_05_15,
                  Silt_15_30, Silt_30_60, ThetaR_00_05, ThetaR_05_15, ThetaR_15_30, ThetaR_30_60, ThetaS_00_05, ThetaS_05_15,
                  ThetaS_15_30, ThetaS_30_60, PH_00_05, PH_05_15, PH_15_30, PH_30_60])
  
// Convert image collection to multi-band image, rename bands
var POLARIS_I = POLARIS_IC.toBands().rename(['BD_00_05', 'BD_05_15', 'BD_15_30', 'BD_30_60', 'Clay_00_05', 'Clay_05_15', 'Clay_15_30', 'Clay_30_60',
                  'OM_00_05', 'OM_05_15', 'OM_15_30', 'OM_30_60', 'Sand_00_05', 'Sand_05_15', 'Sand_15_30', 'Sand_30_60', 'Silt_00_05', 'Silt_05_15',
                  'Silt_15_30', 'Silt_30_60', 'ThetaR_00_05', 'ThetaR_05_15', 'ThetaR_15_30', 'ThetaR_30_60', 'ThetaS_00_05', 'ThetaS_05_15',
                  'ThetaS_15_30', 'ThetaS_30_60', 'PH_00_05', 'PH_05_15', 'PH_15_30', 'PH_30_60'])
  
// Calculate zonal statistics for AIM plots
var fire_pts_polaris = POLARIS_I.reduceRegions({
    collection: fire_pts_all,
    reducer: ee.Reducer.mean(),
    scale: 30
})


// --------------- Export the reduceRegions FeatureCollection to a csv ---------------

Export.table.toDrive({
    collection: fire_pts_polaris,
    description:'Soil_Predictors_Fires',
    fileFormat: 'csv'
});


// ------------ Compile POLARIS images into single image for exporting -------------

var POLARIS_export_I = POLARIS_I.select('PH_30_60', 'OM_30_60')

// Reproject and resample image by nearest neighbors to match Landsat
var LS_ref = ee.Image('users/zackrwerner/landsat_harm_reference')
var LS_proj = LS_ref.projection() // get projection of landsat harmonized image
var POLARIS_I_resample = POLARIS_export_I.reproject(LS_proj) //reproject defaults to nn 

// Calculate focal means at 100 meter scale using reduceNeighborhood
var POLARIS_I_resample_focal = POLARIS_I_resample.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.circle(100, 'meters')})
    .rename(['PH_30_60', 'OM_30_60'])


// Get list of bandnames to export over
var namelist = POLARIS_I_resample_focal.bandNames().getInfo();
print(namelist,'namelist');

for (var i = 0; i < 2; i++) {
// Get image
var img = POLARIS_I_resample_focal.select(namelist[i]).clip(GBbounds.geometry());
    
// Export results to drive
Export.image.toDrive({
        image: img,
        description: namelist[i],
        scale: 30,
        maxPixels: 1e13,
        region: GBbounds});}  