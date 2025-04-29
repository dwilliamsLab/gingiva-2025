//Removing unregistered areas from sample H5 (experimental batch 1 containing 7 Hoechst channels). Not using Channel 18 which is the Hoechst from cycle 5 which has an artifact. Registration is efficient and the rest of the image is not affected

//Creating folders and how to process an entire folder

inputFolder = getDirectory("Input directory");
SaveTarget1 = inputFolder + "Unregistered_areas" + File.separator;
SaveTarget2 = getDirectory("Output directory");
if ( !(File.exists(SaveTarget1)) ) { File.makeDirectory(SaveTarget1); }
processFolder(inputFolder);

 
function processFolder(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if(File.isDirectory(list[i]))
            processFolder("" + input + list[i]);
        if(endsWith(list[i], "t.ims"))
            processFile(input, list[i]);
    }
}
 
function processFile(input, file) {

 //Open output.ims file, make composite, duplicate and split channels
open(inputFolder + file);
run("Make Composite", "display=Composite");
run("Split Channels");
//select overlapping Hoechsts (regions present in all Cycles)
//Briefly, process each Hoechst with background subtraction, median filter and then threshold it. Select the mask and dilate it twice. The reason that I am dilating it, is in the event that there is slight difference that is only observed in the borders of the Hoechsts (nuclei overlapping in more than 80-90%).
//After adding the masks of dilated Hoechsts in ROI manager, use the "AND" function to pick overlapping areas. No you have the ROI that is present in all Hoechst's.
selectImage("C1-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
selectImage("C5-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
selectImage("C9-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
selectImage("C13-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
selectImage("C16-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
selectImage("C23-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Mean dark no-reset");
run("Create Mask");
run("Dilate");
wait(500);
run("Dilate");
wait(500);
run("Dilate");
run("Create Selection");
roiManager("Add");
roiManager("Select", newArray(0,1,2,3,4,5));
roiManager("AND");
roiManager("Add");
roiManager("Select", newArray(0,1,2,3,4,5));
roiManager("Delete");
selectImage("C1-output-1.ims");
close();
selectImage("C5-output-1.ims");
close();
selectImage("C9-output-1.ims");
close();
selectImage("C13-output-1.ims");
close();
selectImage("C16-output-1.ims");
close();
selectImage("C23-output-1.ims");
close();
roiManager("Select", 0);
roiManager("Save", SaveTarget1 + "overlapping.roi");
roiManager("Deselect");
//Find non-overlapping regions from the different Hoechsts 
//This will be done by subtracting overlapping ROI from each of the Hoechst masks. Then what is remaining is the non-overlapping Hoechst. 
//C1
selectImage("C1-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Triangle dark");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C1unregistered");
close();
//C5
selectImage("C5-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
setAutoThreshold("Triangle dark");
run("Create Mask");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C5unregistered");
close();
//C9
selectImage("C9-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
setAutoThreshold("Triangle dark");
run("Create Mask");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C9unregistered");
close();
//C13
selectImage("C13-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
setAutoThreshold("Triangle dark");
run("Create Mask");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C13unregistered");
close();
//C16
selectImage("C16-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
setAutoThreshold("Triangle dark");
run("Create Mask");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C16unregistered");
close();
//C23
selectImage("C23-output.ims");
run("Duplicate...", " ");
run("Median...", "radius=2");
run("Subtract Background...", "rolling=50");
roiManager("Select", 0);
run("Clear", "slice");
roiManager("Show All");
roiManager("Show None");
setAutoThreshold("Triangle dark");
run("Create Mask");
run("Analyze Particles...", "size=120-Infinity pixel show=Masks");
run("Fill Holes");
run("Create Selection");
roiManager("Add");
saveAs("Tiff", SaveTarget1+"C23unregistered");
close();
//Merge non-overlapping Hoechst regions and remove them from all channels.
//After finding all the non-overlapping ROIs from the different cycles, add them with the "Combine" function. This way you have the ROI of all non-overlapping areas, that you can then delete from all channels in the multichannel tif.
//Delete Hoechsts from the first 5 cycles and CD68 channel (poor performance).
//Merge all channels in one multichannel tif file that will be saved in a different folder.
roiManager("Select", newArray(1,2,3,4,5,6));
roiManager("Combine");
roiManager("Add");
roiManager("Select", newArray(0,1,2,3,4,5,6));
roiManager("Delete");
selectImage("C1-output-1.ims");
close();
selectImage("C5-output-1.ims");
close();
selectImage("C9-output-1.ims");
close();
selectImage("C13-output-1.ims");
close();
selectImage("C16-output-1.ims");
close();
selectImage("C23-output-1.ims");
close();
run("Images to Stack", "  title=mask use");
close();
selectImage("C1-output.ims");
close();
selectImage("C5-output.ims");
close();
selectImage("C9-output.ims");
close();
selectImage("C13-output.ims");
close();
selectImage("C16-output.ims");
close();
selectImage("C18-output.ims");
close();
selectImage("C10-output.ims");
close();
run("Images to Stack", "  title=C use");
roiManager("Select", 0);
run("Clear", "stack");
run("Stack to Hyperstack...", "order=xyczt(default) channels=19 slices=1 frames=1 display=Composite");
roiManager("Select", 0);
wait(500);
roiManager("Save", SaveTarget1 + "non-overlapping.roi");
roiManager("Delete");
roiManager("Show All");
roiManager("Show None");
//Manually crop outside the tissue area (based on remaining overlapping Hoechst).
//Even though all unregistered nuclei have already been removed (which is sufficient since segmentation is requiring the presence of nucleus), with this way you crop and keep the "bulk" area that has all the residual tissue and remove major unregistered areas completely (and not only in the nuclei regions).
//This is a way to remove for example a major epithelial area (expressing CK5 from the first cycle) that was later detached. The nuclei had already been removed but the other chennels were not.
////select ROI of tissue and click outside. If multiple ROIs combine them into one and remove them to only keep ONE combined ROI
waitForUser ("select remaining tissue, then hit OK");
roiManager("select", 0);
run("Clear Outside", "stack");
roiManager("Select", 0);
roiManager("Save", SaveTarget1 + "cropped.roi");
roiManager("Delete");
roiManager("Show All");
roiManager("Show None");
//Change order of CD3 and Thy-1 manually (they were in different order in the first experimental batch). Stacks<Tools<Stack Sorter. They are Channel 10 and 11.
waitForUser ("Change order, then hit OK");
//Ideally create a string to name them based on the name of the InputFolder, but I did not find out how to do it, so I just need to manually rename
saveAs("Tiff", SaveTarget2+"H5");
close();






