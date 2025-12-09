//set channels + POI name
POI_name = "Slk19";
channels = newArray("DIC", "KTs", POI_name, "SPBs");

//set channel colors below. If there is no channel for a particular color, put 0.
channel_colors = newArray(

    "KTs",      //red
    0,          //green
    0,          //blue
    "DIC",      //gray
    "SPBs",     //cyan
    0,          //magenta
    POI_name    //yellow
    
);
enhance_contrast_sat = "0.005";

Dialog.create("Choose Directory");
Dialog.addMessage("Choose the directory containing DeltaVision files (channels + Z-stacks) with multiple cells.");
Dialog.show();
dir = getDir("Choose a Directory");

filenames = getFileList(dir); //returns an array with list of files in this folder. We only want TIFFs or .dv files.

saving_folder = 'single_cells'; //where to save the single cell images.
single_cell_macro_datafiles_folder = 'single_cell extra macro files';
RGB_images_folder = 'RGB_images';
ROIs_folder = 'ROIs';

if (File.isDirectory(dir + "/" + saving_folder) == 1) {
    Dialog.create("Are you sure?");
    Dialog.addMessage("It appears you already picked single cells in this folder.\nContinuing will only open images that you have not\npicked cells from yet.");
    Dialog.show();
}

info_file = single_cell_macro_datafiles_folder + "/already_analyzed.txt";
strain_info_file = single_cell_macro_datafiles_folder + "/strain_info.txt";
RGB_images_folder = single_cell_macro_datafiles_folder + "/RGB_images";
ROIs_folder = single_cell_macro_datafiles_folder + "/ROIs";
ran_macro_before = 1; //default: the user has already run this macro.

if (!File.exists(dir + info_file)) {
    //the user has not run the macro before. Let's create the folders for saving the images.
    File.makeDirectory(dir + single_cell_macro_datafiles_folder);
    File.makeDirectory(dir + saving_folder);
    File.makeDirectory(dir + RGB_images_folder);
    File.makeDirectory(dir + ROIs_folder);
    ran_macro_before = 0; //the user has not run this macro.
}

//functions:

//get the index of a value in an array.
function index(value, arr) {
    return_val = -1;
    for (i=0; i < lengthOf(arr); i++) {
        if (arr[i] == value) {
            return_val = i; // Return the index if a match is found
        }
    }
    
    return return_val;
}

//reads the "SBY" strain number from an image title containing underscores (_) as spacers.

function get_strain_name(image_title) {
    return_val = 0;
    strain_name_arr = split(image_title, "_");
    for (z = 0; z < lengthOf(strain_name_arr); z++)  {
        if (startsWith(strain_name_arr[z], "SBY")) {
            return_val = strain_name_arr[z];
            break;
        }
    }
    
    return return_val;
}

//if the user has already run the macro, it returns the number of single cells already exported from a strain.
    //return: 0 if the strain has not been analyzed; >0 = cell count.
function return_cell_count(strain) {
    
    strain_file_string = File.openAsString(dir + "/" + strain_info_file);

    if (strain_file_string != "") {
        
        lines = split(strain_file_string, "\n");

        //parse the tabs.
        strain_name_array = newArray();
        number_of_cells_array = newArray();

        //get the values in the strain_file

        for (line_num = 0; line_num < lengthOf(lines); line_num++) {
            tab_delimiter = split(lines[line_num], "\t");
            strain_saved = tab_delimiter[0];
            cell_count = tab_delimiter[1];

            strain_name_array = Array.concat(strain_name_array, strain_saved);
            number_of_cells_array = Array.concat(number_of_cells_array, cell_count);
        }

        //get the current strain index.
        strain_index = index(strain, strain_name_array);
        
        if (strain_index != -1) {
            
            return_val = round(parseFloat(number_of_cells_array[strain_index])); //integer.
            
        } else {
         
            return_val = 0;
            
        }
        
        return return_val; //this is the number stored in the file
        
    } else {
     
        return 0;
        
    }
}

//adds a 1 to the current cell count in the strain_info file.
function update_cell_count(strain) {
    
    strain_name_array = newArray();
    number_of_cells_array = newArray();
     
    //update the strain_info_file.
    strain_file_string = File.openAsString(dir + "/" + strain_info_file);

    if (strain_file_string != "") {

        lines = split(strain_file_string, "\n");

        for (line_num = 0; line_num < lengthOf(lines); line_num++) {
            tab_delimiter = split(lines[line_num], "\t");
            strain_saved = tab_delimiter[0];
            cell_count = tab_delimiter[1];

            strain_name_array = Array.concat(strain_name_array, strain_saved);
            number_of_cells_array = Array.concat(number_of_cells_array, cell_count);
        }

        //get the current strain index.
        strain_index = index(strain_name, strain_name_array);
        
        if (strain_index != -1) {
           
            current_count = number_of_cells_array[strain_index];
            current_count = parseFloat(current_count); //almost an integer.
            current_count = round(current_count); //now it's an integer.

            number_of_cells_array[strain_index] = current_count + 1;
            
        } else {
        
            //this strain has not been added to this file yet.
            strain_name_array = Array.concat(strain_name_array, strain);
            number_of_cells_array = Array.concat(number_of_cells_array, 1);
            
        }

    } else {
        
        //there is no info so this must be the first entry.
        strain_name_array[0] = strain_name;
        number_of_cells_array[0] = 1;

    }
    
    strain_info_string = "";
    
    for (i = 0; i < lengthOf(strain_name_array); i++) {
        strain_info_string = strain_info_string + strain_name_array[i] + "\t" + number_of_cells_array[i];
        if (i != lengthOf(strain_name_array)-1) {
            strain_info_string = strain_info_string + "\n";    
        }
    }
    
    //lets overwrite the file now with the updated info.
    strain_filehandle = File.open(dir + "/" + strain_info_file);
    print(strain_filehandle, strain_info_string);
    File.close(strain_filehandle);
}

function deleteAllROIs() {
 
    all_rois = newArray();
    
    if (roiManager("Count") != 0) {
     
        for (r = 0; r < roiManager("Count"); r++) {
            all_rois = Array.concat(all_rois, r);            
        }

        roiManager("select", all_rois);
        roiManager("delete"); //delete all ROIs.
        
    }
    
}

function convert_cell_number( cell_num ) {
    
    cell_num_string = "" + cell_num;
    
    if (lengthOf(cell_num_string) == 1) {
        converted_num = "00" + cell_num_string;
    }
    
    if (lengthOf(cell_num_string) == 2) {
        converted_num = "0" + cell_num_string;
    }
    
    if (lengthOf(cell_num_string) == 3) {
        converted_num = cell_num_string;
    }
    
    return converted_num;
    
}

deleteAllROIs();

if (lengthOf(filenames) != 0) {
    
    for (img_ID=0; img_ID < lengthOf(filenames); img_ID++) {
        
        //look to see if there is an information file for this macro (i.e., has the user run this macro in the specified folder before?)
        if (File.exists(dir + "/" + info_file)) {
            //the user has run the macro before. Let's skip all the files we have already analyzed.
            already_analyzed_files = File.openAsString(dir + "/" + info_file);
            files_and_cell_numbers = split(already_analyzed_files, "\n");
            already_analyzed_files_array = newArray();
            for (f=0; f<lengthOf(files_and_cell_numbers); f++) {
                already_analyzed_filenames = split(files_and_cell_numbers[f], '\t');
                already_analyzed_files_array = Array.concat(already_analyzed_files_array, already_analyzed_filenames[0]);
            }
            
            if (index(filenames[img_ID], already_analyzed_files_array) >= 0) {
                //the file has already been analyzed. Skip this file and "continue" to the next.  
                continue;     
            }
        } else {
            //the file does not exist, let's create it.
            info_file_handle = File.open(dir + "/" + info_file);
            File.close(info_file_handle);
            strain_info_filehandle = File.open(dir + "/" + strain_info_file);
            File.close(strain_info_filehandle);
            
        }
        
        setBatchMode(true);
        
        //make sure the file is .dv or .tif:
        check_file = split(filenames[img_ID], ".");
        file_extension = check_file[lengthOf(check_file)-1];
        if (file_extension == "dv") {
            run("Bio-Formats Windowless Importer", "open=[" + dir + "/" + filenames[img_ID] + "]");         
        } else if (file_extension == "tif") {
            open(dir + "/" + filenames[img_ID]);   
        } else {
            continue; //skip the rest of the code and go to the next file.       
        }
        
        //let's duplicate the image, outline the cells via thresholding, save those as ROIs, then make a Z-projection and merge the fluorescence data.
        selectImage(1);
        image_name = getTitle();

        //get the strain name:
        strain_name = get_strain_name(image_name);

        run("Duplicate...", "duplicate");
        selectImage(2);
        rename(image_name + "_selections");

        //this is a hyperstack (i.e. channel AND Z information). Need to use the following function.
        //open the DIC channel, move to the middle Z slice to get a clear representation of cells.
        num_of_slices_per_channel = Math.round((nSlices/lengthOf(channels))/2);
        Stack.setPosition((index("DIC", channels)+1), num_of_slices_per_channel, 1);

        //the background for each cell perimiter is very dark. we can threshold this and make general outlines of the cells to help with visualization.
        /*setThreshold(0, 9670, "raw");
        run("Create Selection");
        run("ROI Manager...");
        roiManager("Add");
        run("Select None");
        setThreshold(65535, 65535, "raw"); //un-threshold stuff so they dont look red.
        */
        //make Z-projection:
        run("Z Project...", "projection=[Max Intensity]");
        selectImage(2); //delete the stack we used to threshold.
        close();
        selectImage(2); //split the channels
        run("Split Channels");

        for(i=1; i<=lengthOf(channels); i++) {
            selectImage(i+1);
            rename(channels[i-1]);
        }

        //selectImage("DIC"); //delete this channel. Don't need, now that we outlined the cells.
        //close();

        for(i=3; i<=5; i++) {
            selectImage(i);
            run("Subtract Background...", "rolling=50"); //subtract background for the user.
            run("Enhance Contrast...", "saturated=" + enhance_contrast_sat);
        }

        merge_channels_string = "";
        possible_channels = newArray('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7');
        for(c=0; c < lengthOf(channel_colors); c++) {

            if (channel_colors[c] != 0) {
                if (merge_channels_string != "") {
                    merge_channels_string = merge_channels_string + " "; //add a space between if it's not the first item.
                }
                merge_channels_string = merge_channels_string + possible_channels[c] + "=" + channel_colors[c];    
            }

        }

        run("Merge Channels...", merge_channels_string);
        //the new name of this image is "RGB";
        //get title without .tif.
        
        title_without_ext = split(image_name, ".");
        RGB_image_name = "";
        if (lengthOf(title_without_ext)>1) {
            for(i=0;i<lengthOf(title_without_ext);i++) {
                if (i != (lengthOf(title_without_ext)-1)) {
                    RGB_image_name = RGB_image_name + title_without_ext[i];  
                }
            }
        } else {
            RGB_image_name = image_name;
        }
        
        RGB_image_name = RGB_image_name + "_RGB.tif";
        selectImage("RGB");
        rename(RGB_image_name);
        saveAs("TIFF", dir + "/" + RGB_images_folder + "/" + RGB_image_name); //save the RGB image.
        
        //next, go to the ROI manager and select the ROI that outlined the cells.
        /*selectImage(RGB_image_name);
        run("Select None");
        roiManager("Set Fill Color", "white");
        roiManager("Select", 0);
        roiManager("Fill");
        roiManager("Select", 0);
        roiManager("Delete");*/
        
        setBatchMode("show");
        
        setTool("rotrect");
        roiManager("Show All with labels");
        
        Dialog.createNonBlocking("Box selections");
        Dialog.addMessage("Outline each cell using the rotating rectangle tool.\nClick \"Add\" in the ROI manager for each cell.\nTry not to overlap cells.\n\nPress OK when done.");
        Dialog.show();
        
        cell_number = return_cell_count(strain_name) + 1;
        starting_cell_number = cell_number;
        
        for (r = 0; r < roiManager("Count"); r++) {
            
            run("Select None");
            roiManager("deselect");
            selectImage(image_name);
            roiManager("select", r);
            roiName = "" + cell_number;
            roiManager("rename", roiName);
            new_image_name = strain_name + "_CELL" + convert_cell_number(cell_number) + "_" + image_name;
            run("Duplicate...", "title=" + new_image_name + " duplicate");
            saveAs("Tiff", dir + "/" + saving_folder + "/" + new_image_name);
            close(); 
            
            update_cell_count(strain_name);
            cell_number = cell_number + 1; //the new number for the next loop.
            
        }
        
        last_cell_num = cell_number-1;
        
        roiManager("Save", dir + "/" + ROIs_folder + "/" + image_name + "_RoiSet__CELLS_" + convert_cell_number(starting_cell_number) + "-" + convert_cell_number(last_cell_num) + ".zip");
        close("*");
        
        deleteAllROIs();
        
        //add this image name to the info_file so we do not analyze it again if the user re-runs the macro in this folder.
        info_file_string = image_name;
        File.append(info_file_string, dir + "/" + info_file);
        
    }
    
    Dialog.create("Macro complete!");
    Dialog.addMessage("All images have been selected and cropped.\nThey are now in the 'single_cells' folder.\nIf you want to restart the macro on the folder,\ndelete all the macro files and folders\nfrom this directory.");
    Dialog.show();
    
} else {

Dialog.create("No images");
Dialog.addMessage("There are no images in this folder.\nPlease run the macro and try again.");
Dialog.show();
    
}