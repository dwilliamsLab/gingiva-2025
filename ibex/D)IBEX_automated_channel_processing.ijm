//Creating folders and how to process an entire folder

inputFolder = getDirectory("Registered_files");
SaveTarget1 = getDirectory("Processed_files");
SaveTarget2 = getDirectory("Masks");
SaveTarget2a = SaveTarget2 + "CK19" + File.separator;
if ( !(File.exists(SaveTarget2a)) ) { File.makeDirectory(SaveTarget2a);
SaveTarget2b = SaveTarget2 + "CK5" + File.separator;
if ( !(File.exists(SaveTarget2b)) ) { File.makeDirectory(SaveTarget2b);
SaveTarget2c = SaveTarget2 + "MCT" + File.separator;
if ( !(File.exists(SaveTarget2c)) ) { File.makeDirectory(SaveTarget2c);
SaveTarget2d = SaveTarget2 + "CD31" + File.separator;
if ( !(File.exists(SaveTarget2d)) ) { File.makeDirectory(SaveTarget2d);
SaveTarget2e = SaveTarget2 + "CD4" + File.separator;
if ( !(File.exists(SaveTarget2e)) ) { File.makeDirectory(SaveTarget2e);
SaveTarget2f = SaveTarget2 + "PanCK" + File.separator;
if ( !(File.exists(SaveTarget2f)) ) { File.makeDirectory(SaveTarget2f);
SaveTarget2g = SaveTarget2 + "MPO" + File.separator;
if ( !(File.exists(SaveTarget2g)) ) { File.makeDirectory(SaveTarget2g);
SaveTarget2h = SaveTarget2 + "Ki67" + File.separator;
if ( !(File.exists(SaveTarget2h)) ) { File.makeDirectory(SaveTarget2h);
SaveTarget2i = SaveTarget2 + "CD8a" + File.separator;
if ( !(File.exists(SaveTarget2i)) ) { File.makeDirectory(SaveTarget2i);
SaveTarget2j = SaveTarget2 + "Thy-1" + File.separator;
if ( !(File.exists(SaveTarget2j)) ) { File.makeDirectory(SaveTarget2j);
SaveTarget2k = SaveTarget2 + "CD3" + File.separator;
if ( !(File.exists(SaveTarget2k)) ) { File.makeDirectory(SaveTarget2k);
SaveTarget2l = SaveTarget2 + "CD45" + File.separator;
if ( !(File.exists(SaveTarget2l)) ) { File.makeDirectory(SaveTarget2l);
SaveTarget2m = SaveTarget2 + "S100a8-9" + File.separator;
if ( !(File.exists(SaveTarget2m)) ) { File.makeDirectory(SaveTarget2m);
SaveTarget2n = SaveTarget2 + "CD138" + File.separator;
if ( !(File.exists(SaveTarget2n)) ) { File.makeDirectory(SaveTarget2n);
SaveTarget2o = SaveTarget2 + "aSMA" + File.separator;
if ( !(File.exists(SaveTarget2o)) ) { File.makeDirectory(SaveTarget2o);
SaveTarget2p = SaveTarget2 + "Hoechst" + File.separator;
if ( !(File.exists(SaveTarget2p)) ) { File.makeDirectory(SaveTarget2p);
SaveTarget2q = SaveTarget2 + "CD20" + File.separator;
if ( !(File.exists(SaveTarget2q)) ) { File.makeDirectory(SaveTarget2q);
SaveTarget2r = SaveTarget2 + "HLA-DR" + File.separator;
if ( !(File.exists(SaveTarget2r)) ) { File.makeDirectory(SaveTarget2r);
SaveTarget2s = SaveTarget2 + "Vimentin" + File.separator;
if ( !(File.exists(SaveTarget2s)) ) { File.makeDirectory(SaveTarget2s);

list = getFileList(inputFolder);
for (i = 0; i < list.length; i++) {
       
 //Open all registered files from the folder, rename them, and split channels
open(inputFolder + list[i]);
   rename("registered");
   run("Split Channels");
   
//Select each individual channel and process accordingly. The median radius selection is based on a regular object size (eg larger for epithelial cell markers) as well as on how "noisy" each channel is (eg slightly bigger in CD20 and CD4). Additionally the BG subtraction rolling ball radius is based on the object size.
//The thresholding method depends on the relative signal/noise ratio and will be selected by an expert histopathologist to include the signal that captures truly positive cells and eliminates background.
//After selecting the thresholded area we generated a mask that included positive pixels (above threshold) but also eliminated small objects (smaller than the minimum truly positive cell, eg larger cut-off for epithelial cell markers)
//The aforementioned mask will be superimposed to the original image and every signal outside of the ROI will be removed. 

//Process CK19
selectImage("C1-registered");
run("Duplicate...", " ");
run("Median...", "radius=1");
run("Subtract Background...", "rolling=100");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=250-Infinity pixel show=Masks");
selectImage("C1-registered-1");
close();
selectImage("Mask of C1-registered-1");
selectImage("Mask of C1-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C1-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C1-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C1-registered-1");
saveAs("Tiff", SaveTarget2a + list[i]);
close();
//Process CK5
selectImage("C2-registered");
run("Duplicate...", " ");
run("Median...", "radius=6");
run("Subtract Background...", "rolling=100");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=50-Infinity pixel show=Masks");
selectImage("C2-registered-1");
close();
selectImage("Mask of C2-registered-1");
selectImage("Mask of C2-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C2-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C2-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C2-registered-1");
saveAs("Tiff", SaveTarget2b + list[i]);
close();
//Process MCT
selectImage("C3-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Moments dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C3-registered-1");
close();
selectImage("Mask of C3-registered-1");
selectImage("Mask of C3-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C3-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C3-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C3-registered-1");
saveAs("Tiff", SaveTarget2c + list[i]);
close();
//Process CD31
selectImage("C4-registered");
run("Duplicate...", " ");
run("Median...", "radius=4");
run("Subtract Background...", "rolling=40");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C4-registered-1");
close();
selectImage("Mask of C4-registered-1");
selectImage("Mask of C4-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C4-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C4-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C4-registered-1");
saveAs("Tiff", SaveTarget2d + list[i]);
close();
//Process CD4
selectImage("C5-registered");
run("Duplicate...", " ");
run("Median...", "radius=4");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C5-registered-1");
close();
selectImage("Mask of C5-registered-1");
selectImage("Mask of C5-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C5-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C5-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C5-registered-1");
saveAs("Tiff", SaveTarget2e + list[i]);
close();
//Process PanCK
selectImage("C6-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=100");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=500-Infinity pixel show=Masks");
selectImage("C6-registered-1");
close();
selectImage("Mask of C6-registered-1");
selectImage("Mask of C6-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C6-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C6-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C6-registered-1");
saveAs("Tiff", SaveTarget2f + list[i]);
close();
//Process MPO
selectImage("C7-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=50-Infinity pixel show=Masks");
selectImage("C7-registered-1");
close();
selectImage("Mask of C7-registered-1");
selectImage("Mask of C7-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C7-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C7-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C7-registered-1");
saveAs("Tiff", SaveTarget2g + list[i]);
close();
//Process Ki67
selectImage("C8-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=50-Infinity pixel show=Masks");
selectImage("C8-registered-1");
close();
selectImage("Mask of C8-registered-1");
selectImage("Mask of C8-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C8-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C8-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C8-registered-1");
saveAs("Tiff", SaveTarget2h + list[i]);
close();
//Process CD8a
selectImage("C9-registered");
run("Duplicate...", " ");
run("Median...", "radius=3");
run("Subtract Background...", "rolling=5");
setAutoThreshold("Moments dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C9-registered-1");
close();
selectImage("Mask of C9-registered-1");
selectImage("Mask of C9-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C9-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C9-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C9-registered-1");
saveAs("Tiff", SaveTarget2i + list[i]);
close();
//Process Thy-1
selectImage("C10-registered");
run("Duplicate...", " ");
run("Median...", "radius=10");
run("Subtract Background...", "rolling=10");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=150-Infinity pixel show=Masks");
selectImage("C10-registered-1");
close();
selectImage("Mask of C10-registered-1");
selectImage("Mask of C10-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C10-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C10-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C10-registered-1");
saveAs("Tiff", SaveTarget2j + list[i]);
close();
//Process CD3
selectImage("C11-registered");
run("Duplicate...", " ");
run("Median...", "radius=3");
run("Subtract Background...", "rolling=20");
setAutoThreshold("Otsu dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C11-registered-1");
close();
selectImage("Mask of C11-registered-1");
selectImage("Mask of C11-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C11-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C11-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C11-registered-1");
saveAs("Tiff", SaveTarget2k + list[i]);
close();
//Process CD45
selectImage("C12-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C12-registered-1");
close();
selectImage("Mask of C12-registered-1");
selectImage("Mask of C12-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C12-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C12-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C12-registered-1");
saveAs("Tiff", SaveTarget2l + list[i]);
close();
//Process S100a8/9
selectImage("C13-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Moments dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C13-registered-1");
close();
selectImage("Mask of C13-registered-1");
selectImage("Mask of C13-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C13-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C13-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C13-registered-1");
saveAs("Tiff", SaveTarget2m + list[i]);
close();
//Process CD138
selectImage("C14-registered");
run("Duplicate...", " ");
run("Median...", "radius=4");
run("Subtract Background...", "rolling=80");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=80-Infinity pixel show=Masks");
selectImage("C14-registered-1");
close();
selectImage("Mask of C14-registered-1");
selectImage("Mask of C14-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C14-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C14-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C14-registered-1");
saveAs("Tiff", SaveTarget2n + list[i]);
close();
//Process aSMA
selectImage("C15-registered");
run("Duplicate...", " ");
run("Median...", "radius=4");
run("Subtract Background...", "rolling=80");
setAutoThreshold("Triangle dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=80-Infinity pixel show=Masks");
selectImage("C15-registered-1");
close();
selectImage("Mask of C15-registered-1");
selectImage("Mask of C15-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C15-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C15-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C15-registered-1");
saveAs("Tiff", SaveTarget2o + list[i]);
close();
//Process Hoechst
selectImage("C16-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=50-Infinity pixel show=Masks");
selectImage("C16-registered-1");
close();
selectImage("Mask of C16-registered-1");
selectImage("Mask of C16-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C16-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C16-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C16-registered-1");
saveAs("Tiff", SaveTarget2p + list[i]);
close();
//Process CD20
selectImage("C17-registered");
run("Duplicate...", " ");
run("Median...", "radius=4");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Moments dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=200-Infinity pixel show=Masks");
selectImage("C17-registered-1");
close();
selectImage("Mask of C17-registered-1");
selectImage("Mask of C17-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C17-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C17-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C17-registered-1");
saveAs("Tiff", SaveTarget2q + list[i]);
close();
//Process HLA-DR
selectImage("C18-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C18-registered-1");
close();
selectImage("Mask of C18-registered-1");
selectImage("Mask of C18-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C18-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C18-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C18-registered-1");
saveAs("Tiff", SaveTarget2r + list[i]);
close();
//Process Vimentin
selectImage("C19-registered");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=60");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");
run("Analyze Particles...", "size=60-Infinity pixel show=Masks");
selectImage("C19-registered-1");
close();
selectImage("Mask of C19-registered-1");
selectImage("Mask of C19-registered-1");
run("Create Selection");
getRawStatistics(nPixels, mean);
print(mean);
if(mean>0){
roiManager("Add");
selectImage("C19-registered");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Select", 0);
roiManager("Delete");
}
else{
selectImage("C19-registered");	
getDimensions(width, height, channels, slices, frames);
print(width, height);
makeRectangle(0, 0, width, height);
run("Clear");
}
selectImage("Mask of C19-registered-1");
saveAs("Tiff", SaveTarget2s + list[i]);
close();
//Merge them in one processed file
run("Images to Stack", "  title=C use");
run("Stack to Hyperstack...", "order=xyczt(default) channels=19 slices=1 frames=1 display=Composite");
saveAs("Tiff", SaveTarget1+ list[i]);
run("Close All");
}


