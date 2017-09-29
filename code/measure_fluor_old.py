import numpy
import scipy.ndimage as ndimage
import freeimage
import pathlib
from concurrent import futures



def hysteresis_threshold(image, low_threshold, high_threshold):
    """Find all regions where image > low_threshold and which contain at least
    one pixel > high_threshold."""
    high_mask = image >= high_threshold
    low_mask = image > low_threshold
    return ndimage.binary_dilation(high_mask, mask=low_mask, iterations=-1)

def fill_small_radius_holes(mask, max_radius):
    outside = ndimage.binary_dilation(numpy.zeros_like(mask), mask=~mask, iterations=-1, border_value=1)
    holes = ~(mask | outside)
    large_hole_centers = ndimage.binary_erosion(holes, iterations=max_radius)
    large_holes = ndimage.binary_dilation(large_hole_centers, mask=holes, iterations=-1)
    small_holes = holes ^ large_holes
    return mask | small_holes

def get_background(mask, offset_radius, background_radius):
    offset = ndimage.binary_dilation(mask, iterations=offset_radius)
    background = ndimage.binary_dilation(offset, iterations=background_radius)
    return background ^ offset

def get_areas(mask):
    labels, num_regions = ndimage.label(mask)
    region_indices = numpy.arange(1, num_regions + 1)
    areas = ndimage.sum(numpy.ones_like(mask), labels=labels, index=region_indices)
    return labels, region_indices, areas

def get_largest_object(mask):
    labels, region_indices, areas = get_areas(mask)
    largest = region_indices[areas.argmax()]
    return labels == largest

def get_mask(image, low_pct, high_pct, max_hole_radius):
    low_thresh, high_thresh = numpy.percentile(image, [low_pct, high_pct])
    mask = hysteresis_threshold(image, low_thresh, high_thresh)
    mask = fill_small_radius_holes(mask, max_hole_radius)
    mask = get_largest_object(mask)
    
    return mask

def measure_intensity(image):
    mask = get_mask(image, 99.2, 99.9, 12)
    background = get_background(mask, 15, 20)
    background_value = numpy.percentile(image[background], 25)
    pixel_data = image[mask] - background_value
    return (mask.sum(), pixel_data.sum())+ tuple(numpy.percentile(pixel_data, [50, 95]))

def measure_intensities(image_files):
    all_values = []
    for image_file in image_files:
        print(str(image_file))
        image = freeimage.read(image_file)
        all_values.append(measure_intensity(image))
    #threadpool = futures.ThreadPoolExecutor(max_workers=4)
    #all_values = numpy.array(list(threadpool.map(measure_intensity, sorted(image_files))))

    return all_values

def get_well_names(image_files):
    wells = []
    for image_file in (image_files):
        row = image_file.name[0] # first letter is row, then ' - ', then two-digit col
        if image_file.name[2]=='_':
            col='0'+image_file.name[1]
        else:    
            col=image_file.name[1:3]
            
        well = row + col
        wells.append(well)
        #wells=sorted(wells)
    return wells

def write_intensities(all_values, wells, csv_out):
    data_header = ['well', 'area', 'integrated', 'median', '95th']
    data = [data_header]
    for well, values in zip(wells, all_values):
        
        data.append(map(str, [well] + list(values)))
    outdata = '\n'.join(','.join(row) for row in data)
    with open(csv_out, 'w') as f:
        f.write(outdata)
    return wells

def holly(image_dir, dye_name, strain_name):
    image_dir = pathlib.Path(image_dir)
    assert image_dir.exists()
    image_files = list(image_dir.glob('*'+dye_name+'*'))
    wells = get_well_names(image_files)
    all_values = measure_intensities(image_files)
    write_intensities(all_values, wells, strain_name+ ' ' + dye_name+'.csv')

def travis(image_dir):
    image_dir = pathlib.Path(image_dir)
    assert image_dir.exists()
    image_files = list(image_dir.glob('*FITC*'))
    wells = get_well_names(image_files)
    image_data = measure_intensities(image_files)
    area = image_data[:,0]
    good_worms = area > 2000
    image_data = image_data[good_worms]
    integrated = image_data[:,1]
    wells = [w for w, g in zip(wells, good_worms) if g]
    num_to_take = int(round(len(wells) * 0.1))
    integ_order = integrated.argsort()
    high_wells = sorted([wells[i] for i in integ_order[-num_to_take:]])
    low_wells = sorted([wells[i] for i in integ_order[:num_to_take]])
    print('High GFP:')
    for w in high_wells:
        print(' '+w)
    print('Low GFP:')
    for w in low_wells:
        print(' '+w)