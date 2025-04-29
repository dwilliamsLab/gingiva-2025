//Creating folders and how to process an entire folder

inputFolder = getDirectory("Registered_files");
SaveTarget1 = getDirectory("Output_directory/segmentation_files");

list = getFileList(inputFolder);
for (i = 0; i < list.length; i++) {
       
 //Open all registered files from the folder, rename them, and split channels
open(inputFolder + list[i]);
   rename("registered");
   run("Split Channels");
// Add the brightest pixels between the first two cell border marker channels. For this multichannel image channel 1 is CK19 and channel 2 is CK5. With this apporach the brightest pixels between those two channels are selected to generate "Result of C1-registered"
imageCalculator("Add create", "C1-registered","C2-registered");
wait(500);
//Pick the brightest pixels between "Result of C1 and the next cell border marker (channel 3- MCT). The resulting image will be "Result of Result of C1".
imageCalculator("Add create", "Result of C1-registered","C3-registered");
wait(500);
//Iterate through all channels with the exception of Ki-67 and Hoechst (nuclei markers.
imageCalculator("Add create", "Result of Result of C1-registered","C4-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of C1-registered","C5-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of C1-registered","C6-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of C1-registered","C7-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of C1-registered","C9-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of C1-registered","C10-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C11-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C12-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C13-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C14-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C15-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C17-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C18-registered");
wait(500);
imageCalculator("Add create", "Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered","C19-registered");
wait(500);
//Merge the nuclear channel (Hoechst in our case) with the final Result of imagecalculator("Result of Result... of C1"). In the subsequent example 17 markers were used for this composite so the word Result is repeated 16 times.
run("Merge Channels...", "c1=C16-registered c2=[Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of Result of C1-registered] create keep");
run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=1 frames=1 display=Grayscale");
saveAs("Tiff", SaveTarget1 + list[i]);
run("Close All");
}