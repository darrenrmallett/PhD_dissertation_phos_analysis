import svg, time, proteins, sys, pprint as pp

export_folder = 'output'
unix_timestamp = int(time.time())
unix_timestamp = str(unix_timestamp)

def generate(all_phosphosites, unix_timestamp_of_xlsx_file) :
    
    protein_sites = {}
    
    for protein, sites in all_phosphosites.items():
        
        #delete the proteins that do not have phosphorylation. 'sites' is a dictionary. It should be empty if there are no phos.
        
        if len(sites) > 0:
            protein_sites[protein] = list(sites.keys())
    
    filehandle = open('output/' + unix_timestamp_of_xlsx_file + '.svg', 'w')

    svg_objects = []

    num = 1
    protein_widths_list = []
    
    #scaling info:
    x_scale = 1.5 #scale = # AAs * 1.5 = pixel width
    start_x = 100
    start_y = 100

    for protein, sites in protein_sites.items():

        protein_info = proteins.Protein(protein)

        x_width = protein_info.length*x_scale
        protein_widths_list.append(x_width)

        y_height = 20 #always 20 px)
        start_y = start_y+(y_height*3) #3-times the y position, times by the num.

        x_text_start = start_x-60
        y_text_start = start_y+15

        #protein name:
        svg_objects.append( svg.Text(x=x_text_start, y=y_text_start, style="font: 15px sans-serif;", text=protein) )

        #protein rectangle object:
        svg_objects.append( svg.Rect(x=start_x, y=start_y, width=x_width, height=y_height, fill="#B9E3FF", stroke="#006FB0", stroke_width=1) )

        #phosphosites:

        for site in sites:

            x_position = (site * x_scale) + start_x
            display_text = protein_info.subsequence(site, site+1) + str(site)

            svg_objects.append( svg.Line(x1 = x_position, x2 = x_position, y1 = start_y-10, y2 = start_y, stroke="#000000", stroke_width=1.5) )            
            svg_objects.append( svg.Text(x = x_position-5, y = start_y-15, style="font: 8px sans-serif;", text=display_text) )            

        num+=1

        final_height = start_y

    max_width = max(protein_widths_list)

    if max_width > 2200:

        canvas_width = max_width + 200
    else:
        canvas_width = 2200

    canvas = svg.SVG(
        width=canvas_width,
        height=final_height+100,
        elements = svg_objects
    )

    canvas = str(canvas)

    filehandle.write(canvas)

    filehandle.close()