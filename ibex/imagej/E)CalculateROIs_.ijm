inputFolder = getDirectory("Registered_files");
SaveTarget1 = getDirectory("Area_calculation");
SaveTarget2 = getDirectory("Area_coordinates");
SaveTarget3 = getDirectory("Area_masks");

list = getFileList(inputFolder);

for (fileIndex = 0; fileIndex < list.length; fileIndex++) {
    // Open all registered files from the folder, rename them, and split channels
    open(inputFolder + list[fileIndex]);
    rename("registered");
    run("Split Channels");

    // Process Vim
    selectImage("C19-registered");
    run("Duplicate...", " ");
    run("Median...", "radius=5");
    run("Subtract Background...", "rolling=5");
    setAutoThreshold("Percentile dark no-reset");
    run("Analyze Particles...", "size=200-Infinity pixel show=Masks");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Create Selection");
    roiManager("add");

    // Process CK5
    selectImage("C2-registered");
    run("Duplicate...", " ");
    run("Median...", "radius=5");
    run("Subtract Background...", "rolling=10");
    setAutoThreshold("Percentile dark no-reset");
    run("Analyze Particles...", "size=500-Infinity pixel show=Masks");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Create Selection");
    roiManager("add");

    // Process PanCK
    selectImage("C6-registered");
    run("Duplicate...", " ");
    run("Median...", "radius=5");
    run("Subtract Background...", "rolling=10");
    setAutoThreshold("Percentile dark no-reset");
    run("Analyze Particles...", "size=500-Infinity pixel show=Masks");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Create Selection");
    roiManager("add");

    // Process CD138
    selectImage("C14-registered");
    run("Duplicate...", " ");
    run("Median...", "radius=5");
    run("Subtract Background...", "rolling=10");
    setAutoThreshold("Percentile dark no-reset");
    run("Analyze Particles...", "size=500-Infinity pixel show=Masks");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Dilate");
    run("Create Selection");
    roiManager("add");

    // Merge Epi ROIs
    roiManager("Select", newArray(1, 2, 3));
    roiManager("Combine");
    roiManager("add");
    roiManager("Select", newArray(1, 2, 3));
    roiManager("delete");

    // Find Epi ROI
    roiManager("Select", newArray(0, 1));
    roiManager("AND");
    roiManager("add");
    roiManager("Select", newArray(1, 2));
    roiManager("XOR");
    roiManager("add");
    roiManager("Select", newArray(1, 2));
    roiManager("delete");

    selectImage("C16-registered");
    run("Duplicate...", " ");
    setAutoThreshold("Percentile dark");
    roiManager("Select", 1);
    run("Clear Outside");
    run("Create Mask");
    run("Analyze Particles...", "size=3000-Infinity pixel show=Masks");
    run("Create Selection");
    roiManager("Add");
    roiManager("Select", 1);
    roiManager("delete");

    // Find tissue ROI
    roiManager("Select", newArray(0, 1));
    roiManager("combine");
    roiManager("add");
    // Find Crevicular and Surface ROI, if no crevicular select don't select anything
    waitForUser("select Area containing surface and add to manager, then hit OK");
    roiManager("deselect");
    nROIs = roiManager("count");

    // Use an if statement to perform different actions based on the number of ROIs
    if (nROIs < 4) {
        // Save Masks of Epi and CT
        roiManager("Select", 0);
        roiManager("Set Fill Color", "cyan");
        roiManager("Select", 1);
        roiManager("Set Fill Color", "red");
        selectImage("Mask of C19-registered-1");
        roiManager("Show All");
        run("Flatten");
        saveAs("Tiff", SaveTarget3 + list[fileIndex] + "EPI_CT");

        // Find immune ROI
        selectImage("C12-registered");
        run("Duplicate...", " ");
        run("Median...", "radius=5");
        run("Subtract Background...", "rolling=20");
        setAutoThreshold("Triangle dark no-reset");
        run("Analyze Particles...", "size=100-Infinity pixel show=Masks");
        run("Create Selection");
        roiManager("add");
        roiManager("Select", newArray(0, 3));
        roiManager("AND");
        roiManager("add");
        roiManager("deselect");
        roiManager("Select", 3);
        roiManager("delete");

        // Measure area and save measurements
        roiManager("Select", 0);
        run("Measure");
        roiManager("Select", 1);
        run("Measure");
        roiManager("Select", 2);
        run("Measure");
        roiManager("Select", 1);
        run("Measure");
        setResult("Area", nResults(), 0);
        roiManager("Select", 3);
        run("Measure");
        setResult("Category", 0, "CT");
        setResult("Category", 1, "Epi");
        setResult("Category", 2, "Total");
        setResult("Category", 3, "Surface");
        setResult("Category", 4, "Crevicular");
        setResult("Category", 5, "Immune");
        selectWindow("Results");
        saveAs("txt", SaveTarget1 + list[fileIndex] + ".csv");
        selectWindow("Results");
        run("Clear Results");

        // Save coordinates
        roiManager("Select", 0);
        roiManager("Rename", "CT");
        roiManager("Select", 1);
        roiManager("Rename", "Epi/Surface");
        roiManager("Select", 2);
        roiManager("Rename", "Total");
        roiManager("Select", 3);
        roiManager("Rename", "Immune");

        numROIs = roiManager("count");
        nr = 0;

        // Define mapping between numerical labels and desired labels
        labels = newArray("CT", "Epi/Surface", "Total", "Immune");

        for (roiIndex = 0; roiIndex < numROIs; roiIndex++) {
            roiManager("Select", roiIndex);
            Roi.getCoordinates(x, y);

            for (coordIndex = 0; coordIndex < x.length; coordIndex++) {
                setResult("Label", coordIndex + nr, labels[roiIndex]); // Set the label to the desired label
                setResult("X", coordIndex + nr, x[coordIndex]);
                setResult("Y", coordIndex + nr, y[coordIndex]);
            }
            nr += x.length;
            updateResults();
        }
        saveAs("txt", SaveTarget2 + list[fileIndex] + ".csv");
        run("Clear Results");

        roiManager("Select", newArray(0, 1, 2, 3));
        roiManager("save", SaveTarget2 + list[fileIndex]+ ".zip");
        roiManager("deselect");

        // Save masks
        roiManager("Select", newArray(1, 2));
        roiManager("delete");
        selectImage("Mask of C19-registered-1");
        roiManager("Select", 1);
        roiManager("Set Fill Color", "magenta");
        roiManager("Show All");
        run("Flatten");
        saveAs("Tiff", SaveTarget3 + list[fileIndex] + "immuneROI");
        run("Close All");
        roiManager("Select", newArray(0, 1));
        roiManager("delete");
        run("Clear Results");
    } else {
        roiManager("Select", newArray(1, 3));
        roiManager("AND");
        roiManager("add");
        roiManager("deselect");
        roiManager("Select", 3);
        roiManager("delete");
        roiManager("Select", newArray(1, 3));
        roiManager("XOR");
        roiManager("add");

        // Save Masks of Epi and CT
        roiManager("Select", 0);
        roiManager("Set Fill Color", "cyan");
        roiManager("Select", 3);
        roiManager("Set Fill Color", "red");
        roiManager("Select", 4);
        roiManager("Set Fill Color", "orange");
        selectImage("Mask of C19-registered-1");
        roiManager("Show All");
        run("Flatten");
        saveAs("Tiff", SaveTarget3 + list[fileIndex] + "EPI_CT");

        // Find immune ROI
        selectImage("C12-registered");
        run("Duplicate...", " ");
        run("Median...", "radius=5");
        run("Subtract Background...", "rolling=20");
        setAutoThreshold("Triangle dark no-reset");
        run("Analyze Particles...", "size=100-Infinity pixel show=Masks");
        run("Create Selection");
        roiManager("add");
        roiManager("Select", newArray(0, 5));
        roiManager("AND");
        roiManager("add");
        roiManager("deselect");
        roiManager("Select", 5);
        roiManager("delete");

        // Measure area and save measurements
        roiManager("Select", 0);
        run("Measure");
        roiManager("Select", 1);
        run("Measure");
        roiManager("Select", 2);
        run("Measure");
        roiManager("Select", 3);
        run("Measure");
        roiManager("Select", 4);
        run("Measure");
        roiManager("Select", 5);
        run("Measure");
        setResult("Category", 0, "CT");
        setResult("Category", 1, "Epi");
        setResult("Category", 2, "Total");
        setResult("Category", 3, "Surface");
        setResult("Category", 4, "Crevicular");
        setResult("Category", 5, "Immune");
        selectWindow("Results");
        saveAs("txt", SaveTarget1 + list[fileIndex] + ".csv");
        selectWindow("Results");
        run("Clear Results");

        // Save coordinates
        roiManager("Select", 0);
        roiManager("Rename", "CT");
        roiManager("Select", 1);
        roiManager("Rename", "Epi");
        roiManager("Select", 2);
        roiManager("Rename", "Total");
        roiManager("Select", 3);
        roiManager("Rename", "Surface");
        roiManager("Select", 4);
        roiManager("Rename", "Crevicular");
        roiManager("Select", 5);
        roiManager("Rename", "Immune");

        numROIs = roiManager("count");
        nr = 0;

        // Define mapping between numerical labels and desired labels
        labels = newArray("CT", "Epi", "Total", "Surface", "Crevicular", "Immune");

        for (roiIndex = 0; roiIndex < numROIs; roiIndex++) {
            roiManager("Select", roiIndex);
            Roi.getCoordinates(x, y);

            for (coordIndex = 0; coordIndex < x.length; coordIndex++) {
                setResult("Label", coordIndex + nr, labels[roiIndex]); // Set the label to the desired label
                setResult("X", coordIndex + nr, x[coordIndex]);
                setResult("Y", coordIndex + nr, y[coordIndex]);
            }
            nr += x.length;
            updateResults();
        }
        saveAs("txt", SaveTarget2 + list[fileIndex] + ".csv");
        run("Clear Results");

        roiManager("Select", newArray(0, 1, 2, 3, 4, 5));
        roiManager("save", SaveTarget2 + list[fileIndex]+ ".zip");
        roiManager("deselect");

        // Save masks
        roiManager("Select", newArray(1, 2, 3, 4));
        roiManager("delete");
        selectImage("Mask of C19-registered-1");
        roiManager("Select", 1);
        roiManager("Set Fill Color", "magenta");
        roiManager("Show All");
        run("Flatten");
        saveAs("Tiff", SaveTarget3 + list[fileIndex] + "immuneROI");
        run("Close All");
        roiManager("Select", newArray(0, 1));
        roiManager("delete");
        run("Clear Results");
    }
}
