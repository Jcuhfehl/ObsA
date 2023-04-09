from astropy.wcs import WCS
import astropy.io.ascii as asc
import math
import astroalign as aa
import astropy.io.fits as fits
import ccdproc
import matplotlib.pyplot as plt
import numpy as np
import imexam
import os


def table(rows):
    for column in rows:
        for element in column:
            print("{}".format(element).rjust(20), end="")
        print()


table_data = [["Name", "Filter", "Exposure time", "Mean counts"]]
for j in range(1, 11):
    filename = '../flats/flat%dV.FIT' % j
    data, header = fits.getdata(filename, header=True)
    row = [filename, header["FILTER"], header["EXPTIME"], np.mean(data)]
    table_data.append(row)

table(table_data)


def add_pedestal(filenames, output_filenames):
    for j in range(len(filenames)):
        file_in = filenames[j]
        file_out = output_filenames[j]

    # - Read the header and image data and extract the PEDESTAL card

        data = fits.getdata(file_in)
        hdr = fits.getheader(file_in)
        pedestal = 1.0*hdr['PEDESTAL']

        print("Adding %d ADU: %s -> %s" % (pedestal, file_in, file_out))

    # Add back the PEDESTAL

        datap = data + pedestal
        hdr['PEDESTAL'] = 0  # Update the PEDESTAL keyword for completeness' sake
        hdr['COMMENT'] = 'PEDESTAL of %d ADU added' % pedestal

    # Save the result

        fits.writeto(file_out, datap, hdr, overwrite=True)


def inv_avg(data):
    # Function that returns the inverse average of the values
    # in the data array
    return 1/data.mean()


def create_masterflat(filenames, output_filename, skip_pedestal=False):
    if not skip_pedestal:
        temp_filenames = [os.path.basename(filename) for filename in filenames]
        add_pedestal(filenames, temp_filenames)

    # Make a list of the filenames
    masterflat = ccdproc.combine(temp_filenames, method='average', unit='adu',
                                 scale=inv_avg)
    masterflat.write(output_filename, overwrite=True)
    return masterflat


flatV = create_masterflat(['../flats/flat%dV.FIT' %
                          j for j in range(1, 11)], 'flatV.fits')
flatr = create_masterflat(['../flats/flat%dr.FIT' %
                          j for j in range(1, 11)], 'flatr.fits')
# Calculate mean pixel value in each combined flat-field,
# should be equal to 1


print("Mean level of combined V masterflat = {}".format(np.mean(flatV)))
print("Mean level of combined r masterflat = {}".format(np.mean(flatr)))

sdv = np.std(flatV)
sdr = np.std(flatr)

print(sdv, sdr)


# flatV.write("flatV.fits")
# flatr.write("flatr.fits")


def flatfield_correction(filenames_in, filenames_out, masterflat_filename):
    # For each science image:
    for j in range(len(filenames_in)):
        file_in = filenames_in[j]
        file_out = filenames_out[j]
        print("Processing file: %s -> %s" % (file_in, file_out))
        data, header = fits.getdata(file_in, header=True)
        masterflat = fits.getdata(masterflat_filename)

        # - Add back the PEDESTAL to the image data and divide by the masterflat
        pedestal = header["PEDESTAL"]
        datap = data + pedestal
        dataf = datap / masterflat

        # - Add comments to the header
        # Update the PEDESTAL keyword for completeness' sake
        header['PEDESTAL'] = 0
        header['COMMENT'] = 'PEDESTAL of %d ADU added and divided by the flatfield' % pedestal

        # Save the flat-fielded science image together with the updated header
        fits.writeto(file_out, dataf, header, overwrite=True)


filenames_in = ["../science/science%dV.FIT" % j for j in range(1, 11)]
filenames_out = ["science%dVf.fits" % j for j in range(1, 11)]
flatfield_correction(filenames_in, filenames_out, "flatV.fits")

filenames_in = ["../science/science%dr.FIT" % j for j in range(1, 11)]
filenames_out = ["science%drf.fits" % j for j in range(1, 11)]
flatfield_correction(filenames_in, filenames_out, "flatr.fits")


def align_images(source_filename, target_filename, output_filename):
    datatarget = fits.getdata(target_filename)
    datasrc = fits.getdata(source_filename)
    T, (source_pos_array, target_pos_array) = aa.find_transform(datasrc, datatarget)
    print("Image name = {}  shifts= {} rot= {}  scale = {}"
          .format(source_filename, T.translation, T.rotation*180/math.pi, T.scale))
    # Ensure the rotation is not too big
    assert (abs(T.rotation*180/math.pi) < 0.01)
    assert (abs(T.scale - 1) < 0.0002)  # Ensure the scale is not to different
    data_tran, footprint = aa.apply_transform(
        T, np.int32(datasrc), np.int32(datatarget))

    fits.writeto(output_filename, data_tran, overwrite=True)

# Align the flat-fielded science images:


def align_multiple_images(target_filename, source_filenames, output_filenames):
    for i in range(len(source_filenames)):
        source_filename = source_filenames[i]
        output_filename = output_filenames[i]
        align_images(source_filename, target_filename, output_filename)


def combine_and_align_images(science_filenames, align_target, output_filename):
    temp_filenames = ["temp{}.fits".format(
        i) for i in range(len(science_filenames))]
    align_multiple_images(align_target, science_filenames, temp_filenames)
    sciVavg = ccdproc.combine(temp_filenames, method="average", unit="adu")
    print(np.array(sciVavg).mean())

    # Save the output
    sciVavg.write('tmpcmb.fits', overwrite=True)
    # ...
    cmb, hdr = fits.getdata('tmpcmb.fits', 0, header=True)
    fits.writeto(output_filename, cmb, hdr, overwrite=True)
    print('Done!')


align_target = "science1Vf.fits"
filenames_V = ["science{}Vf.fits".format(i) for i in range(1, 11)]
combine_and_align_images(filenames_V, align_target, "Vfinal.fits")

filenames_R = ["science{}rf.fits".format(i) for i in range(1, 11)]
combine_and_align_images(filenames_R, align_target, "rfinal.fits")


dataV = fits.getdata("Vfinal.fits")
dataR = fits.getdata("rfinal.fits")

mean_V = np.mean(dataV)
std_V = np.std(dataV)
mean_R = np.mean(dataR)
std_R = np.std(dataR)


os.system("sex Vfinal.fits -CATALOG_NAME V.cat")
os.system("sex Vfinal.fits,rfinal.fits -CATALOG_NAME r.cat")


# Read the data

Vdata = asc.read('V.cat')
rdata = asc.read('r.cat')


name, ra, de = [], [], []
with open("urat.txt") as f:
    for line in f:
        lsplit = line.split()
        name.append(lsplit[0])
        ra.append(float(lsplit[1]))
        de.append(float(lsplit[2]))

w = WCS("new-image.fits")
xy = w.all_world2pix(ra, de, 0)
print(w)

pcoo = []
for _x, _y, _name in zip(xy[0], xy[1], name):
    pcoo.append((_x, _y, _name))

# Read the data
Vdata = asc.read('V.cat')
rdata = asc.read('r.cat')


def match_star_magnitudes(urat_name):
    """This function matches a star in the catalogue generated by␣
    ,→sourceextractor to the urat catalogue by comparing the distance between the␣
    ,→coordinates between the two"""
    pcoo_element = None
    for element in pcoo:
        pcoo_element = element
        if element[2] == urat_name:
            star_urat = element
            break
    star_Vdata = None
    for star in Vdata:
        if np.linalg.norm(np.array(star[1], star[2]) - np.array(star_urat[0], star_urat[1])) < 3:
            star_Vdata = star
            star_rdata = rdata[star_Vdata[0]-1]
            break
    if star_Vdata is None:
        return []
    return urat_name, star_Vdata[5], star_Vdata[6], star_rdata[5], star_rdata[6], pcoo_element


urat = []
urat_with_instrumental_magnitudes = []
pcoo_new = []
with open("urat.txt") as f:
    for line in f:
        lsplit = line.split()
        urat.append(lsplit)
for star in urat:
    matched = match_star_magnitudes(star[0])
    if matched == [] or matched[2] > 0.01 or matched[4] > 0.01:
        continue
    if len(star) != 7:
        print("Oh no")
        continue
    new_star = star
    for data in matched[1:]:
        new_star.append(str(data))
    print(f"Matched {new_star[0]}")
    urat_with_instrumental_magnitudes.append(new_star)
    pcoo_new.append(matched[-1])

V_offsets = []
r_offsets = []
for star in urat_with_instrumental_magnitudes:
    V_offsets.append(float(star[3])-float(star[7]))
    r_offsets.append(float(star[5])-float(star[9]))
V_offsets = np.array(V_offsets)
r_offsets = np.array(r_offsets)
V_offset = np.mean(V_offsets)
r_offset = np.mean(r_offsets)
print("V offset:", V_offset)
print("r offset:", r_offset)


Vdata = asc.read('V.cat')
rdata = asc.read('r.cat')

for i in range(len(Vdata)):  # Add the zero point value to the source extractor data
    Vdata[i]["MAG_APER"] += V_offset
    rdata[i]["MAG_APER"] += r_offset

# Now make the plot. Set the x- and y limits, add labels to the axes
g_r = []
g = []
for i in range(len(Vdata)):
    if Vdata[i]["MAGERR_APER"] > 1 or rdata[i]["MAGERR_APER"] > 1:
        continue  # I had quite a lot of really inaccurate data, this gets rid of it
    if Vdata[i]["X_IMAGE"] > 2000 and Vdata[i]["Y_IMAGE"] < 220:
        continue  # This is the area of the field that had ice on it
    g_r.append(Vdata[i]["MAG_APER"]-rdata[i]["MAG_APER"])
    g.append(Vdata[i]["MAG_APER"])
g_r = np.array(g_r)
g = np.array(g)

B_V = 0.90*g_r + 0.21  # https://www.sdss4.org/dr12/algorithms/sdssubvritransform/
V = g - 0.58*g_r - 0.01
V_absolute_offset = -9.9  # This is manually determined


MVH, BVH = np.loadtxt('Hyades.txt', usecols=(1, 2), unpack=True)
MVP, BVP = np.loadtxt('Pleiades.txt', usecols=(1, 2), unpack=True)

# isochrones = np.loadtxt("isochrones.dat", usecols=(2, 28, 29))
isochrones = np.loadtxt("isochrones2.dat", usecols=(0, 8, 9))
isochrones_split = y = [isochrones[isochrones[:, 0] == k]
                        for k in np.unique(isochrones[:, 0])]


def apparent_magnitude_cmd():
    plt.plot(g_r, g, '.', label="M36")  # This is in the g and r filters
    plt.legend()
    plt.xlim(-0.3, 1.8)
    plt.ylim(np.max(g)+1, np.min(g)-1)
    plt.xlabel(r'$g-r$')
    plt.ylabel(r'$M_g$')
    plt.title("CMD of M36 with apparent magnitudes")
    plt.minorticks_on()
    plt.show()


def absolute_magnitude_cmd():
    # This is in the g and r filters
    plt.plot(g_r, g+V_absolute_offset, '.', label="M36")
    plt.legend()
    plt.xlim(-0.3, 1.8)
    plt.ylim(12, -4)
    plt.xlabel(r'$g-r$')
    plt.ylabel(r'$M_g$')
    plt.title("CMD of M36 with absolute magnitudes")
    plt.minorticks_on()
    plt.show()


def absolute_magnitude_cmd_compare():
    plt.plot(BVH, MVH, '.', label='Hyades')
    plt.plot(BVP, MVP, '.', label='Pleiades')
    plt.plot(B_V, V+V_absolute_offset, '.', label="M36")
    plt.legend()
    plt.xlim(-0.3, 1.8)
    plt.ylim(12, -4)
    plt.xlabel(r'$B-V$')
    plt.ylabel(r'$M_V$')
    plt.title("CMD of M36 with absolute magnitudes")
    plt.minorticks_on()
    plt.show()


def absolute_magnitude_cmd_isochrones():
    #for isochrone in isochrones_split[::10]:
    isochrone = isochrones_split[42]
    xs = isochrone.transpose()[1]-isochrone.transpose()[2]
    ys = isochrone.transpose()[1]
    plt.plot(
        xs, ys, label=f"{10**isochrone.transpose()[0][0]/1000000} Myr")

    plt.plot(B_V, V+V_absolute_offset, '.', label="M36")
    # plt.plot(BVP, MVP, '.', label="Pleiades")
    # plt.plot(BVH, MVH, '.', label="Hyades")

    plt.legend()
    plt.xlim(-0.3, 1.8)
    plt.ylim(12, -4)
    plt.xlabel(r'$B-V$')
    plt.ylabel(r'$M_V$')
    plt.title("CMD of M36 with absolute magnitudes")
    plt.minorticks_on()
    plt.show()


def errors_histograms():
    g_errors = [star["MAGERR_APER"]
                for star in Vdata if star["MAGERR_APER"] < 1]
    r_errors = [star["MAGERR_APER"]
                for star in rdata if star["MAGERR_APER"] < 1]
    g_discarded = len(Vdata)-len(g_errors)
    r_discarded = len(rdata)-len(r_errors)
    print(
        f"The number of discarded stars in g and r respectively: {g_discarded}, {r_discarded}")
    plt.hist([g_errors, r_errors], bins=25, label=["g", "r"])
    plt.legend()
    plt.show()


def show_imexam():
    os.environ["XPA_METHOD"] = "local"
    viewer = imexam.connect()
    new_image = fits.open("Vfinal.fits")
    viewer.load_fits(new_image)
    viewer.scale()
    viewer.mark_region_from_array(pcoo_new, size=20, textoff=25)
    viewer.imexam()


def compare_check_imexam():
    os.environ["XPA_METHOD"] = "local"
    viewer = imexam.connect()
    final = fits.open("rfinal.fits")
    check = fits.open("check.fits")

    viewer.frame(1)
    viewer.load_fits(final)
    viewer.scale(90)
    viewer.zoomtofit()

    viewer.frame(2)
    viewer.load_fits(check)
    viewer.scale(90)
    viewer.zoomtofit()
    viewer.blink(blink=True, interval=10)
    viewer.imexam()


# compare_check_imexam()
# show_imexam()
apparent_magnitude_cmd()
absolute_magnitude_cmd()
absolute_magnitude_cmd_compare()
absolute_magnitude_cmd_isochrones()
# errors_histograms()

D = 10 * 10**(-V_absolute_offset/5)
MSTO = -0.40  # This is in the g filter
age = 10**((MSTO+27.75)/3.18)
age2 = (10**10)*(10**((MSTO-5.23)/-2.5))**(-5/7)
median_counts_in_empty_area_g = [
    49.22640681892971, 45.58912374706513, 45.212821536115655,
    45.544330316146976, 45.665797635717475, 45.04771921972416,
    48.44699409183238, 47.40807889647495, 45.89445556625974, 48.3720535961707]
median_counts_in_empty_area_r = [
    61.30271598596723, 61.836990635782456, 61.04314572837647,
    62.066029894320124, 61.04956179784504, 61.76169194236403,
    61.018569477387175, 62.07226163778414, 61.360477362077745,
    61.91087141932828]
sky_brightness_g = -2.5 * \
    np.log10(np.mean(median_counts_in_empty_area_g)) + V_offset
sky_brightness_r = -2.5 * \
    np.log10(np.mean(median_counts_in_empty_area_r)) + r_offset

fwhm_g = [9.64, 10.33, 11.11, 11.59, 11.00]
fwhm_r = [25.53, 17.19, 27.47, 14.45, 19.42]

offset_max = -10.5
offset_min = -9
distance_min = 10*10**(-offset_min/5)
distance_max = 10*10**(-offset_max/5)

age_max = 10**((MSTO-V_absolute_offset+offset_min+27.75)/3.18)
age_min = 10**((MSTO-V_absolute_offset+offset_max+27.75)/3.18)

print(f"FWHM in g and r filters: {np.mean(fwhm_g)} and {np.mean(fwhm_r)}")

print(
    f"The distance to M36 is determined at {D}pc (+{distance_max-D}/-{D-distance_min}), and its age is {(age/1000000)} (+{(age_max-age)/1000000}/-{(age-age_min)/1000000}) Myr.")
print(
    f"The sky brightness in the g and r filter are {sky_brightness_g} and {sky_brightness_r}")
