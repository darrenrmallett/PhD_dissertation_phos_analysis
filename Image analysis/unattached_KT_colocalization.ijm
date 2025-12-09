//Open image stack, select cells to analyze, save each cell as separate image TIF stack.

setFont("SansSerif", 16);

//set your POI name:
POI_name = 'Stu1';
z_stack = true;

//channel settings, in ORDER:
channels = newArray("DIC", "KTs", POI_name, "SPBs");

//manual thresholding (minimum thresholding value): 
    // set to 'STRAIN' for setting the threshold ONCE per strain (cell #1 will be used, less biased but can be less accurate, depending on how consistence cells are).
    // set to 'ALL' for setting the threshold for every image (this will add accuracy at the cost of some bias)
    // set to 'DETERMINE' to have the program determine the lower threshold value (use 1 pixel value above the MAX pixel value of the background calculated per channel) - No bias but prone to lots of measurement error.
    // set to 'SET' and put in threshold information in the 'set_thresholds' array below, in order of each channel in the hyperstack.
    //  FOR DIC image, we do not want to threshold anything, so this is set to 65535 (nothing will be thresholded unless it's pure white!)
    manual_thresholding = 'ALL';
    set_thresholds = newArray(65535, 333, 240, 900);

//python script pathname for final summary data:
python_script = "/Users/dmallett/Dropbox/Graduate School/Biggins Lab - Main/ImageJ Macros/summarize_KT_localization_results.py";

//set scale:
px_per_micron = 15.3846; //this value I got from the calibrated DeltaVision.

//other settings:
min_particle_size = 8; //in square pixels
results_directory = "Results";
images_directory = results_directory + "/Images";
montage_directory = images_directory + "/Montages";
ROIs_directory = images_directory + "/ROIs";
results_file = results_directory + "/results.txt";

setBatchMode(true);
manual_thresholding = manual_thresholding.toUpperCase;

function contains( value, array ) {
    for (i=0; i<lengthOf(array); i++) {
        if ( array[i] == value ) return true;
    }
    return false;
}

function remove_file_extension(title) {
    
    title = split(title, ".");
    
    if (lengthOf(title) >= 2) {
        new_title_array = Array.deleteIndex(title, (lengthOf(title)-1));
    }
    
    new_title = "";
    
    for (i=0; i<lengthOf(new_title_array); i++) {
        new_title = new_title + new_title_array[i];    
    }
    
    return new_title;
    
}

function update_results(dir, results_array) {
 
    //results array consists of the following (no more, and no less!):
    //dir = directory we are working in.
    /*results = newArray(

        strain_name,
        cell_id,
        image_name,
        label,
        roi,
        status,
        area_px,
        area_um,
        mean,
        mean_minus_bg,
        mode,
        median,
        min,
        max,
        bg_mean,
        RawIntDen,
        notes

    );*/

    result_string = "";

    for (result_item = 0; result_item < lengthOf(results_array); result_item++) {

        result_string = result_string + results_array[result_item];

        if (result_item != (lengthOf(results_array)-1)) {
            result_string = result_string + "\t";   
        }

    }

    File.append(result_string, dir + "/" + results_file);
    
}

function get_strain_name(image_title) {
    //the strain name should be the first item in the image title.
    strain_name_arr = split(image_title, "_");
    return strain_name_arr[0];
}

function get_cell_number(image_title) {
    //the cell number should be the second item in the image title.
    strain_name_arr = split(image_title, "_");
    cell_num_arr = split(strain_name_arr[1], "ELL");
    return cell_num_arr[1]; //second item in the new array is the actual number.
}

function check_image_title_format(image_title) {
    //make sure the image_title has SBY#_CELL#_{whatever} format. For whatever reason, this regex does not work using \d (digit), so we are forced to do an integer character set [0-9]
    if (matches(image_title, 'SBY[0-9]+_CELL[0-9]+_.*\.tif')) {
        return 1;  
    } else {
        return 0;
    }
}

function fix_pathnames_for_python(pathname) {
    new_path = replace(pathname, " ", "\\\ "); //replaces spaces with "\ "
    return new_path;
}

function getMinThreshold() {
    
    set_thresholds = newArray();
    
    //we want to open the channel images and have the user select the threshold.
    for (i = 0; i < lengthOf(channels); i++) {
        
        //we don't want to threshold DIC or POI images.
        if (channels[i] != "DIC" && channels[i] != POI_name) {
            
            selectImage(channels[i]);
            setBatchMode("show");
            run("In [+]");
            run("In [+]"); //zoom in twice... so it's easy to see the image.
            setThreshold(800, 65535, "raw"); //random value to get the slider in the middle somewhere.
            run("Threshold...");
            waitForUser(channels[i] + ' threshold: CELL: ' + cell_id, 'Set the lower threshold that properly\nhighlights the ' + channels[i] + '.\nClick OK when done.');
            getThreshold(lower, upper);
            
            min_threshold = lower;
            
        } else {
         
            min_threshold = 65535;
            
        }
        
        set_thresholds = Array.concat(set_thresholds, min_threshold);
        setBatchMode('hide');
        
    }
    
    return set_thresholds;
    
}

//START MACRO here:
close("*"); //close all images if any are open.
Dialog.create("Choose Directory");
Dialog.addMessage("Choose the directory containing images of single cells in Z-stack form.\nNOTE: these filenames must start with SBY#_CELL#_{whatever} or\nthis macro will NOT work!");
Dialog.show();
dir = getDir("Choose a Directory");

filenames = getFileList(dir); //returns an array with list of files in this folder. We only want TIFFs.

if (File.isDirectory(dir + "/Results") == 1) {
    Dialog.create("Are you sure?");
    Dialog.addMessage("It appears this folder has already been analyzed.\nAnalyzing again will overwrite existing analysis. Continue?");
    Dialog.show();
}

if (lengthOf(filenames) != 0) {
    //create a "results" folder, and a folder to save Images and ROIs, and Montages.
    File.makeDirectory(dir + "/" + results_directory);
    File.makeDirectory(dir + "/" + images_directory);
    File.makeDirectory(dir + "/" + montage_directory);
    File.makeDirectory(dir + "/" + ROIs_directory);
    
    //generate results file and print the headers:    
    results_filehandle = File.open(dir + results_file); //creates the results file
    header_elements = newArray(
    
        'genotype',
        'cell_id',
        'image_name',
        'channel',
        'ROI',
        'status',
        'area_px',
        'area_um',
        'mean',
        'mean_minus_bg',
        'mode',
        'median',
        'min',
        'max',
        'bg_mean',
        'RawIntDen',
        'threshold_info',
        'notes'
        
    );
    //generate header string:
    header_string = "";
    for (h = 0; h < lengthOf(header_elements); h++) {
        header_string = header_string + header_elements[h];
        if (h == (lengthOf(header_elements)-1)) {
            header_string = header_string + "\n";
        } else {
            header_string = header_string + "\t";
        }
    }
    print(results_filehandle, header_string);
    File.close(results_filehandle); //close the file.
    
}

strain_and_cell_id_tracker = ''; //will set to strain_name + "_" + cell_id. Used to determine if it's a new cell or not.
strain_tracker = ''; //will set to strain. Used to determine if it's a new strain or not for thresholding purposes.


for (file = 0; file < lengthOf(filenames); file++) {
    
    filename = filenames[file];
    //we only want TIFFs...
    arr = split(filename, ".");
    last_item = lengthOf(arr)-1;
    file_extension = arr[last_item]; 
    
    if (file_extension == "tif") {
        
        errors = newArray(); //errors array. put errors here if they occur.
        
        open(dir + "/" + filename);
        
        image_name      = getTitle();
        image_width     = getWidth();
        image_height    = getHeight();
        
        if (check_image_title_format(image_name) == 0) {
            exit("Make sure all files within the folder\nhave the format: SBY#_CELL#_{whatever}.\nThe macro will now quit."); 
        }
        
        strain_name = get_strain_name(image_name);
        cell_id = get_cell_number(image_name);

        //if Z-stack, then we want to make max-IP.
        if(z_stack == true) {
            run("Z Project...", "projection=[Max Intensity]");      
        }
        
        selectImage(2);
        run("Duplicate...", "duplicate");
        selectImage(3);
        max_IP_image_title = "" + remove_file_extension(image_name) + "_MAX_IP (MEASURED IMAGE)";
        saveAs("tiff", dir + images_directory + "/" + max_IP_image_title);
        close();

        selectImage(1);
        new_image_title = "" + remove_file_extension(image_name) + '_ORIGINAL';
        saveAs("tiff", dir + images_directory + "/" + new_image_title);
        close();

        //rename image = rename(new_name);

        run("Set Scale...", "distance=1 known=1 unit=pixel"); //RESET the scale because I want the measurements in pixels. Then we can change it back later.

        //Set measurements.
        run("Set Measurements...", "area mean modal min integrated median display add redirect=None decimal=5");

        //split the channels to do the thresholding.
        run("Split Channels");

        //Rename each image with what channel it is.
        for(i = 0; i < 4; i++) {

            selectImage(i+1);
            rename(channels[i]);

        }
        
        if (manual_thresholding == 'STRAIN' && strain_tracker != strain_name) {
            
            set_thresholds = getMinThreshold();
            strain_tracker = strain_name;
            
        }
        
        if (manual_thresholding == "ALL") {
            
            set_thresholds = getMinThreshold();
            
        }
        
        //make RGB image.
        run("Merge Channels...", "c1=KTs c5=SPBs c7=" + POI_name + " keep"); //the image is named 'RGB'
        selectImage("RGB");
        rename("Merge");
        image_name_no_ext = remove_file_extension(image_name);
        merge_name = images_directory + "/" + remove_file_extension(image_name) + "_MERGE";
        saveAs("tiff", dir + "/" + merge_name);
        rename("Merge");
        
        //next, we want to get bg fluorescence to (1) measure bg fluorescence, and (2) determine minimum threshold value for particle picking IF the user has it set to 'DETERMINE'
        //we will make a 25px x 25px box in the upper lefthand corner. this is arbitrary, ASSUMING there are no ROIs within this region.
        
        selectImage(channels[0]); //go to DIC channel.
        makeRectangle(0, 0, 25, 25);
        roiManager("Add");
        roiManager("select", 0);
        new_name = "BG";
        bg_ROI_id = 0;
        roiManager("rename", new_name);

        for (i = 0; i < lengthOf(channels); i++) {

            last_ROI = roiManager("count")-1;
            selectImage(channels[i]);
            roiManager("select", last_ROI);
            run("Measure");
            setResult("ROI", (nResults-1), "BG");

        }
        
        if (manual_thresholding == 'SET') {
            
            set_thresholds = newArray();
            
            for (i = 0; i < lengthOf(channels); i++) {
                if (getResultString("ROI", i) == "BG") {
                    confirm_channel = split(getResultLabel(i), ':');
                    if (confirm_channel[0] == channels[i]) {
                        if (confirm_channel[0] == "DIC") {

                            set_thresholds = Array.concat(set_thresholds, 65535);

                        } else {

                            set_thresholds = Array.concat(set_thresholds, round(getResult("Max", i)));

                        }
                    }
                }
            }

        }
        
        //generate the threshold_string for the output file that is used later.
        threshold_string = "";
        for (i = 0; i < lengthOf(set_thresholds); i++) {
            
            if (channels[i] != "DIC" && channels[i] != POI_name) {
             
                threshold_string = threshold_string + channels[i] + " = (" + set_thresholds[i] + ", 65535); ";
                
            }
            
        }
        
        last_result_index = nResults-1;

        //now threshold each image according to threshold settings determined above.

        for(i = 0; i < 4; i++) {

            if (channels[i] != POI_name) {
                //we only want to threshold the SPBs and KTs, and then from there, measure the fluorescence of the POI at those spots.
                selectImage(channels[i]);
                setOption("BlackBackground", true);                
                setThreshold(set_thresholds[i], 65535, "raw");
                roiManager("deselect");
                run("Select None");
                run("Analyze Particles...", "size=" + min_particle_size + "-Infinity pixel show=Overlay display include add composite");
                resetThreshold(); //remove the threshold so it's no longer red.
            }

        }

        //now we have our ROIs. deselect all, if any are selected.
        run("Select None");
        roiManager("deselect");

        //select the ROIs and rename them according to the results table... because the ROI names they are given are funky.
        //Also, put all indices of ROIs for each ROI species into a new array.
        KT_ROIs = newArray(); //the ROI index for KT particles
        SPBs_ROIs = newArray(); //the ROI index for SPB particles
        POI_ROIs = newArray(); //the ROI index for POI particles
        for(i = 1; i < roiManager("count"); i++) {
            //we start at ROI index #1 (the second ROI) since we already have the BG ROI as index 0.

            roiManager("select", i);
            label = getResultLabel(last_result_index + i);

            if (label == "KTs") {
                KT_ROIs = Array.concat(KT_ROIs, i);
                roiManager("rename", "KT_" + i);
                setResult("ROI", last_result_index + i, "KT_" + i);
            }

            if (label == "SPBs") {
                SPBs_ROIs = Array.concat(SPBs_ROIs, i);
                roiManager("rename", "SPBs");
                setResult("ROI", last_result_index + i, "SPBs");
                setResult("status", last_result_index + i, "aKT"); //this is always the case
            }

            if (label == POI_name) {
                POI_ROIs = Array.concat(POI_ROIs, i);
                roiManager("rename", "" + POI_name + (i+1));
                setResult("ROI", last_result_index + i, "" + POI_name + (i+1));
            }

        }

        //make sure there is only ONE SPB focus. More than one indicates this cell still has a spindle, and 
        // hence, the nocodazole treatment did not work - and we don't want to quantify that.
        if (lengthOf(SPBs_ROIs) == 1) {

            //now that they have been renamed according to the results table, we will determine which KT foci are attached or unattached.
                //if they are unattached, they should NOT be associated with the SPBs
                //and, you guessed it... if they are attached, they are associated with the SPBs.

            //A list of ROI indices that are identified as UNATTACHED kinetochores:
            unattached_KT_ROIs = newArray();
            //A list of ROI indices that are identified as ATTACHED kinetochores:
            attached_KT_ROIs = newArray();

            //for each kinetochore, identify whether it co-localizes with SPBs or not:
            for (i = 0; i < lengthOf(KT_ROIs); i++) {

                roiManager("Deselect");
                select_ROIs = Array.concat(KT_ROIs[i], SPBs_ROIs);
                roiManager("Select", select_ROIs);
                roiManager("AND"); //select the area where the two ROIs overlap. No area is selected if they do not overlap.

                //determine if there is an area selected (i.e. did they overlap?)
                if (Roi.size > 0) {
                    //this kinetochore focus is attached to the SPB!            
                    attached_KT_ROIs = Array.concat(attached_KT_ROIs, KT_ROIs[i]);
                    //annotate in the results table that this KT is attached.
                    setResult("status", last_result_index + KT_ROIs[i], "aKT");
                } else {
                    //there is not area selected so, this kinetochore focus is unattached!
                    unattached_KT_ROIs = Array.concat(unattached_KT_ROIs, KT_ROIs[i]);
                    setResult("status", last_result_index + KT_ROIs[i], "uaKT");
                }

            }
            
            //throw an error if there are no unattached kinetochores.
            if (lengthOf(unattached_KT_ROIs) == 0 || lengthOf(attached_KT_ROIs) == 0) {
                
                //there are either (1) no unattached KTs, or (2) no attached KTs. Throw an error for each one.
                
                if (lengthOf(unattached_KT_ROIs) == 0) {
                 
                    errors = Array.concat(errors, "No unattached KTs present.");
                    
                }
                
                if (lengthOf(attached_KT_ROIs) == 0) {
                    
                    errors = Array.concat(errors, "No attached KTs present.");
                    
                }
            
            } else {

                //next, we want to measure the fluorescence intensity of the POI at the attached KTs/SPBs if there is no thresholded signal there for quantification purposes.
                //combine both the SPB and attached KT ROIs.
                roiManager("select", Array.concat(attached_KT_ROIs, SPBs_ROIs));
                roiManager("Combine");
                roiManager("Add");
                roiManager("select", roiManager("count")-1);
                roiManager("rename", "SPB+att_KTs");

                //now measure the POI fluorescence.
                run("Select None");
                selectImage(POI_name);
                roiManager("select", roiManager("count")-1); 
                run("Measure");
                setResult("ROI", (nResults-1), "SPB+att_KTs");
                setResult("status", (nResults-1), "aKT");

                //next, we want to measure the fluorescence intensity of the POI at the unattached KTs.
                for (i=0; i < lengthOf(unattached_KT_ROIs); i++) {

                    selectImage(POI_name);
                    run("Select None");
                    roiManager("select", unattached_KT_ROIs[i]);
                    run("Measure");
                    setResult("ROI", (nResults-1), "KT_" + (unattached_KT_ROIs[i]));
                    setResult("status", (nResults-1), "uaKT");

                }        


                //finally, we want to measure the 'diffuse' signal surrounding the KTs and uaKTs.
                //we need to select a thresholded area that is, ideally, greater than the maximum background fluorescence signal,
                // as measured in our small square. Then, we will EXCLUDE the POI ROIs determined from the initial threshold, so we are not measuring those twice.

                //determine the min threshold we need for this... again, this is equal to the max bg fluorescence signal in our small square.
                for (n = 0; n < nResults; n++) {

                    label = POI_name + ':BG';
                    if (getResultLabel(n) == label) {

                        max_bg_value = getResult('Max', n);
                        set_lower_threshold_val = round(max_bg_value + 1);

                    }

                }
                

                //next, set a threshold with minimum threshold value of max_bg_value + 1
                selectImage(POI_name);
                run("Select None");
                selectImage(POI_name);
                setThreshold(set_lower_threshold_val, 65535, "raw");
                setOption("BlackBackground", true);
                num_of_results_so_far = nResults;
                last_added_ROI = roiManager("count")-1;
                run("Analyze Particles...", "size=200-Infinity pixel show=Overlay display include add composite"); //200 square pixels is good for a minimum of "diffuse" signal. otherwise it will pick up all kinds of small noise particles.
                
                if (roiManager("count")-1 == last_added_ROI) {
                    //there was no ROI added... I.e. we can't measure the diffuse signal for some reason. Throw an error.
                    errors = Array.concat(errors, "Could not measure the diffuse signal.");
                } else {
                
                    //this could have added more than one 'particle'. let's figure out how many it actually added.
                    total_particles = nResults - num_of_results_so_far;

                    //delete these measurements from the results table, we don't want them. we only wanted to make a temporary ROI.
                    IJ.deleteRows(num_of_results_so_far, nResults-1);

                    //combine the 1 or more particles into a single ROI.
                    last_newly_added_ROI = roiManager("count")-1;
                    rois_to_select = newArray();

                    for (l = last_added_ROI+1; l <= last_newly_added_ROI; l++) {
                        rois_to_select = Array.concat(rois_to_select, l);    
                    }

                    run("Select None");
                    roiManager("Deselect");
                    roiManager("select", rois_to_select);
                    roiManager("Combine");
                    roiManager("Add");
                    //now let's delete all those other single particles since we now have one particle combined.
                    roiManager("Deselect");
                    roiManager("select", rois_to_select);
                    roiManager("Delete");                

                    newly_added_ROI = roiManager("count")-1;
                    resetThreshold();

                    //now we want to exclude KT signals, SPB signals, and already quantified POI signals.
                    roiManager("Select", Array.concat(KT_ROIs, SPBs_ROIs, POI_ROIs));
                    roiManager("Combine"); //combine all the signals we already measured into one ROI. Then we can exclude this from the diffuse POI ROI.
                    roiManager("Add");

                    combined_ROI = roiManager("count")-1;
                    roiManager("Select", Array.concat(newly_added_ROI, combined_ROI));
                    roiManager("Combine"); //now this the ENTIRE area of all the quantified ROIs.
                    roiManager("Add");
                    //next, exlusive OR (XOR) the two previous ROIs.

                    roiManager("Select", Array.concat(roiManager("count")-1, roiManager("count")-2));
                    roiManager("XOR"); //we only want the region that is not the union.
                    roiManager("Add"); //now, this is the ROI to quantify diffuse signal.

                    roiManager("select", roiManager("count")-1);
                    roiManager("rename", POI_name + "_diffuse");

                    //delete the temporary ROIs we made.
                    roiManager("select", roiManager("count")-2);
                    roiManager("Delete");
                    roiManager("select", roiManager("count")-2);
                    roiManager("Delete");
                    roiManager("select", roiManager("count")-2);
                    roiManager("Delete");

                    roiManager("select", roiManager("count")-1);
                    roiManager("Update");
                    roiManager("Measure");
                    setResult("ROI", (nResults-1), POI_name + "_diffuse");
                    setResult("status", (nResults-1), "diffuse");


                    //clean up the 'status' column
                    for (i=0; i<nResults; i++) {

                        if (getResult("status", i) == 0) {

                            setResult("status", i, "");

                        }

                        //also set up the 'notes' column.
                        setResult("notes", i, "");

                    }

                    SPBs_bg_r     = "";
                    KTs_bg_r      = "";
                    POI_bg_r      = "";
                    POI_diffuse_r = "";

                    //fetch the results from the results table.
                    for (r = 0; r < nResults; r++) {
                        label = getResultLabel(r);

                        //first, we want to get the background signal for each channel. fetch which rows in the results table these are.
                        is_it_bg = split(label, ':');
                        if (lengthOf(is_it_bg) == 2) {

                            if (is_it_bg[1] == 'BG') {
                                //this is a background measurement.
                                if (is_it_bg[0] == "SPBs") {
                                    SPBs_bg_r = r;
                                }

                                if (is_it_bg[0] == "KTs") {
                                    KTs_bg_r = r;
                                }

                                if (is_it_bg[0] == POI_name) {
                                    POI_bg_r = r;
                                }
                            }

                        }

                        if (lengthOf(is_it_bg) == 2) {

                            if (is_it_bg[1] == 'diffuse') {
                                //this is the POI diffuse measurement. get the row.
                                POI_diffuse_r = r;
                            }

                        }

                    }

                    special_measurements = newArray(SPBs_bg_r, KTs_bg_r, POI_bg_r, POI_diffuse_r);

                    //all the other rows are full measurements we want.
                    //fetch the results from the results table, again.
                    has_colon = newArray();
                    for (r = 0; r < nResults; r++) {
                        label = getResultLabel(r);

                        if (contains(r, special_measurements) == false && startsWith(label, "DIC") == false) {

                            has_colon = split(label, ':');
                            if (lengthOf(has_colon) == 2) {
                                label = has_colon[0];                        
                            }

                            roi     = getResultString("ROI", r);
                            status  = getResultString("status", r);
                            area_px = getResult("status", r);
                            status  = getResultString("status", r);
                            area_px = getResult("Area", r);
                            area_um = sqrt(area_px)/px_per_micron;
                            area_um = pow(area_um, 2);
                            mean    = getResult("Mean", r);
                            mode    = getResult("Mode", r);
                            median  = getResult("Median", r);
                            min     = getResult("Min", r);
                            max     = getResult("Max", r);
                            notes   = getResultString("notes", r);

                            if (label == "SPBs") {
                                bg_mean = getResult("Mean", SPBs_bg_r);
                            }

                            if (label == "KTs") {
                                bg_mean = getResult("Mean", KTs_bg_r);
                            }

                            if (label == POI_name) {
                                bg_mean = getResult("Mean", POI_bg_r);
                            }

                            mean_minus_bg = mean-bg_mean;
                            RawIntDen = getResult("RawIntDen", r);

                            results = newArray(

                                strain_name,
                                cell_id,
                                image_name,
                                label,
                                roi,
                                status,
                                area_px,
                                area_um,
                                mean,
                                mean_minus_bg,
                                mode,
                                median,
                                min,
                                max,
                                bg_mean,
                                RawIntDen,
                                threshold_string,
                                notes

                            );

                            update_results(dir, results);

                        }

                    }


                    //make a montage of all the ROIs.
                    //close('*'); //closes all images.

                    montage_array = newArray("BG_roi", "aKTs_roi", "uaKTs_roi", "SPBs_roi", POI_name + "_diffuse");

                    for (m=0; m < lengthOf(montage_array); m++) {

                        newImage(montage_array[m], "16-bit black", image_width, image_height, 1);

                        run("Select None");
                        roiManager("Deselect");
                        selectImage(montage_array[m]);

                        if (montage_array[m] == "BG_roi") {

                            roiManager("select", bg_ROI_id);
                            setForegroundColor(201, 139, 255);
                            roiManager("Set Fill Color", "purple");
                            roiManager("select", bg_ROI_id);

                        } else if (montage_array[m] == "aKTs_roi") {

                            roiManager("select", attached_KT_ROIs);
                            setForegroundColor(201, 139, 255);
                            roiManager("Set Fill Color", "purple");
                            roiManager("select", attached_KT_ROIs);

                        } else if (montage_array[m] == "uaKTs_roi") {

                            roiManager("select", unattached_KT_ROIs);
                            setForegroundColor(201, 139, 255);
                            roiManager("Set Fill Color", "red");
                            roiManager("select", unattached_KT_ROIs);

                        } else if (montage_array[m] == "SPBs_roi") {

                            roiManager("select", SPBs_ROIs);
                            setForegroundColor(201, 139, 255);
                            roiManager("Set Fill Color", "blue");
                            roiManager("select", SPBs_ROIs);

                        } else if (montage_array[m] == POI_name + "_diffuse") {

                            RoiManager.selectByName(POI_name + "_diffuse");
                            setForegroundColor(201, 139, 255);
                            roiManager("Set Fill Color", "yellow");
                            RoiManager.selectByName(POI_name + "_diffuse");

                        }

                        roiManager("Fill");

                    }

                    run("Images to Stack", "use keep");
                    selectImage("Stack");
                    stack_name = "" + remove_file_extension(image_name) + '_STACK_WITH_ROIs';
                    saveAs("tiff", dir + montage_directory + "/" + stack_name);
                    setForegroundColor(185, 252, 255);
                    run("Make Montage...", "columns=5 rows=2 scale=1 border=2 label use");
                    montage_name = "" + remove_file_extension(image_name) + '_MONTAGE';
                    saveAs("tiff", dir + montage_directory + "/" + montage_name);          

                }
                
            }


        } else {

            //there are more than one SPB signal so there is an issue.

            errors = Array.concat(errors, "More than one SPB detected."); //add an error message for this cell

        }
        
        //select all ROIs and save them.
        //first, get all ROIs into an array.

        ROI_array = newArray();
        for (roi_ID = 0; roi_ID < roiManager("count"); roi_ID++) {
            ROI_array = Array.concat(ROI_array, roi_ID);
        }

        roiManager("select", ROI_array);
        roiManager("Save", dir + ROIs_directory + "/" + remove_file_extension(image_name) + "_ROIs.zip");
        
        roiManager("select", ROI_array);
        roiManager("Delete");

        close('*'); //closes all images.
        run("Clear Results");  //clears the results window.

        if (lengthOf(errors) > 0) {

            //we want to put errors in the 'notes' section.

            label           = "";
            roi             = "";
            status          = "";
            area_px         = "";
            status          = "";
            area_px         = "";
            area_um         = "";
            mean            = "";
            mode            = "";
            median          = "";
            min             = "";
            max             = "";
            notes           = "";
            bg_mean         = "";
            mean_minus_bg   = "";
            RawIntDen       = "";
            notes           = "";

            //concat the errors into one string:
            for (e = 0; e < lengthOf(errors); e++) {
                notes = notes + errors[e];
                if (e != lengthOf(errors)-1) {
                    notes = notes + "; ";   
                }
            }

            results = newArray(

                strain_name,
                cell_id,
                image_name,
                label,
                roi,
                status,
                area_px,
                area_um,
                mean,
                mean_minus_bg,
                mode,
                median,
                min,
                max,
                bg_mean,
                RawIntDen,
                threshold_string,
                notes

            );

            update_results(dir, results);

        }
        
    }
    
}

//this python file will calculate the summarized_results.txt file using the results file generated above    
    //first argument: the python script file location; second arg: location of results folder.    
exec("python3", python_script, dir + results_directory, POI_name);

Dialog.create("Complete!");
Dialog.addMessage("Macro is now complete.");
Dialog.show();

