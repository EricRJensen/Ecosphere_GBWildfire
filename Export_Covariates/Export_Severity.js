//  Import fire points
var fire_pts = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/Fire_points");

// Import and produce times burned image
var Once = ee.FeatureCollection('users/ericjensen41_default/Thesis/Chapter2/FirePolygons/Fires_BurnedOnce')
            .map(function(f){
            return f.set('TimesBurn', 1)})
var MTOnce = ee.FeatureCollection('users/ericjensen41_default/Thesis/Chapter2/FirePolygons/Fires_BurnedMTOnce')
            .map(function(f){
            return f.set('Year', f.get('Fire_Year'))})
var fires = MTOnce.merge(Once)
var nBurn = fires.reduceToImage(['TimesBurn'], ee.Reducer.mode()).rename('TimesBurn')

// Get list of fires names to iterate over
var fires = fire_pts.aggregate_array('Fire_Name').distinct()

// Import dNBR and RdNBR images
var d_1995 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/1995dNBR').set('Year', 1995)
var d_1999 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/1999dNBR').set('Year', 1999)
var d_2000 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/2000dNBR').set('Year', 2000)
var d_2001 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/2001dNBR').set('Year', 2001)
var d_IC = ee.ImageCollection([d_1995, d_1999, d_2000, d_2001])

var r_1995 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/1995RdNBR').set('Year', 1995)
var r_1999 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/1999RdNBR').set('Year', 1999)
var r_2000 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/2000RdNBR').set('Year', 2000)
var r_2001 = ee.Image('users/ericjensen41_default/Thesis/Covariates/MTBS/2001RdNBR').set('Year', 2001)
var r_IC = ee.ImageCollection([r_1995, r_1999, r_2000, r_2001])

// Extract severity values for each fire
var extract_sev = function(fire){
    var fire_plots = fire_pts.filterMetadata('Fire_Name', 'equals', fire)
    var fire_year = ee.Number(fire_plots.first().get('Year'))
    var d_map = ee.Image(d_IC.filterMetadata('Year', 'equals', fire_year).first()).rename(['dnbr'])
    var r_map = ee.Image(r_IC.filterMetadata('Year', 'equals', fire_year).first()).rename(['rdnbr'])
    var sev_map = d_map.addBands(r_map).addBands(nBurn)

    var sev_rr = sev_map.reduceRegions({
        collection: fire_plots,
        reducer: ee.Reducer.mean(),
        scale:30})

    return(sev_rr)
}

// Run reduceRegions for each fire
var fires0 = extract_sev(fires.get(0))
var fires1 = extract_sev(fires.get(1))
var fires2 = extract_sev(fires.get(2))
var fires3 = extract_sev(fires.get(3))
var fires4 = extract_sev(fires.get(4))
var fires5 = extract_sev(fires.get(5))
var fires6 = extract_sev(fires.get(6))
var fires_rr = fires0.merge(fires1).merge(fires2).merge(fires3).merge(fires4)
                    .merge(fires5).merge(fires6)
print(fires_rr.size())
print(fires_rr.first())

// Export CSV
Export.table.toDrive({
    collection: fires_rr,
    description:'Severity_fires',
    fileFormat: 'CSV'
});