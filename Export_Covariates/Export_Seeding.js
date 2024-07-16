// Import seeding and fire points datasets
var Aerial = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/TreatmentPolygons/LTDL_AerialSeed"),
    Drill = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/TreatmentPolygons/LTDL_DrillSeed"),
    fire_pts = ee.FeatureCollection("users/ericjensen41_default/Thesis/Chapter2/Fire_pointsAll");

// Create images for aerial and drill to reduce regions over
var AerialImg = Aerial.map(function(f){
    return(f.set('Aerial', 1))
}).reduceToImage({
    properties: ['Aerial'],
    reducer: ee.Reducer.first()})

var DrillImg = Drill.map(function(f){
    return(f.set('Drill', 1))
}).reduceToImage({
    properties: ['Drill'],
    reducer: ee.Reducer.first()})

// Combine images
var SeedingImg = AerialImg.addBands(DrillImg).rename(['Aerial', 'Drill'])

// Reduce regions over the new images
var seeding_rr = SeedingImg.reduceRegions({collection: fire_pts, reducer: ee.Reducer.mean(), scale:30})

  // Export CSV
Export.table.toDrive({
    collection: seeding_rr,
    description:'Seeding_fires',
    fileFormat: 'CSV'
});

Map.addLayer(Aerial, {}, 'AerialVector')
Map.addLayer(AerialImg, {}, 'AerialVectorImg')
Map.addLayer(Drill, {}, 'DrillVector')
Map.addLayer(DrillImg, {}, 'DrillVectorImg')