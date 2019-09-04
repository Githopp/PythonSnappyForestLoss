#Sentinel Processing
# Polarisation 
polarisation = 'VH'
# Area of interest - Borneo
wkt = 'POLYGON((116.03248117558246122 -0.58024577195271831, 116.0402848949578356 -1.00792931754335768, 115.51509743598057867 -0.99859593410496927, 115.52202811919585201 -0.56881162348975067, 116.03248117558246122 -0.58024577195271831))'

import warnings
warnings.filterwarnings("ignore")
import os
import sys
import glob
#sys.path.append('/opt/anaconda/bin/')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from snappy import jpy
from snappy import ProductIO
from snappy import GPF
from snappy import HashMap

import gc
import zipfile

def zip_extract():
    data_path = os.getcwd() #get working dir
    #print(os.listdir(data_path))
    global ziplist 
    ziplist = list()
    for file in os.listdir(data_path):
        if file.endswith(".zip"):
            ziplist.append(os.path.join(data_path, file))
    
    for i in ziplist:
        dir = os.path.isdir(os.path.splitext(i)[0] + '.SAFE')
        if dir:
            print("{} already extracted".format(i))

        else:
            with zipfile.ZipFile(i,"r") as izip_elem:
                izip_elem.extractall(data_path)
                print(i,"extracted")
    return ziplist

# Feed with ziplist
def process(ziplist):
    product_set=[]

    for f in ziplist:
        f = os.path.splitext(f)[0] #pretty dirty works only with zip, because filenam is cut.
        f = _process(f)
        
        product_set.append(f)
    
    print("Creating Stack with",len(product_set),"bands")
    if len(product_set) <= 1:
        ##RuntimeError: org.esa.snap.core.gpf.OperatorException: Please select at least two source products
        parameters = HashMap()
        parameters.put('resamplingType', None)
        parameters.put('initialOffsetMethod', 'Orbit')
        parameters.put('extent', 'Master')
        create_stack = GPF.createProduct('CreateStack', parameters, product_set[0])

        #write the stack
        ProductIO.writeProduct(create_stack, 'create_stack.dim', 'BEAM-DIMAP')
        print('Single Product created')
    else:
        parameters = HashMap()
        parameters.put('resamplingType', None)
        parameters.put('initialOffsetMethod', 'Product Geolocation')
        parameters.put('extent', 'Master')
        ## Linear to dB
        create_stack = GPF.createProduct('CreateStack', parameters, product_set)

        #write the stack
        ProductIO.writeProduct(create_stack, 'create_stack.dim', 'BEAM-DIMAP')
        #ProductIO.writeProduct(terrain_correction, s1_identifier + '.tif', "GeoTIFF-BigTIFF")
        #ProductIO.writeProduct(terrain_correction, s1_identifier + '.dim', 'BEAM-DIMAP')
        print('Stack created')

def _process(s1_identifier):
    data_path = os.getcwd()
    #s1_identifier = 'S1A_IW_GRDH_1SDV_20160416T215858_20160416T215923_010852_0103DC_FAC8'
    s1meta = "manifest.safe"
    s1prd = os.path.join(data_path, s1_identifier, s1_identifier + '.SAFE', s1meta)
    #print(s1prd)

    reader = ProductIO.getProductReader("SENTINEL-1")
    product = reader.readProductNodes(s1prd, None)

    print("Bands:\n",list(product.getBandNames()),"Length:",len(product.getBandNames()))


    # ### ThermalNoiseRemoval step

    # HashMap = jpy.get_type('java.util.HashMap')    
    # GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    # 
    # parameters = HashMap()
    # 
    # parameters.put('selectedPolarisations', polarisation)
    # parameters.put('removeThermalNoise', 'true')
    # parameters.put('reIntroduceThermalNoise', 'false')
    # 
    # thermal_noise_removal = GPF.createProduct('ThermalNoiseRemoval', parameters, product)

    # ## Apply-Orbit-File

    # parameters = HashMap()
    # 
    # parameters.put('orbitType', 'Sentinel Precise (Auto Download)')
    # parameters.put('polyDegree', '3')
    # parameters.put('continueOnFail', 'false')
    # 
    # apply_orbit_file = GPF.createProduct('Apply-Orbit-File', parameters, thermal_noise_removal)

    # # Subset
    parameters = HashMap()

    parameters.put('sourceBands', 'Amplitude_{}'.format(polarisation))
    #parameters.put('region', '')
    parameters.put('geoRegion', wkt)
    parameters.put('subSamplingX', '1')
    parameters.put('subSamplingY', '1')
    parameters.put('fullSwath', 'false')
    parameters.put('tiePointGridNames', '')
    parameters.put('copyMetadata', 'true')

    subset = GPF.createProduct('Subset', parameters, product)
    print('Subset {} done'.format(wkt))
    ## Calibration
    parameters = HashMap()

    parameters.put('auxFile', 'Product Auxiliary File')
    parameters.put('outputImageInComplex', 'false')
    parameters.put('outputImageScaleInDb', 'false')
    parameters.put('createGammaBand', 'true')
    parameters.put('createBetaBand', 'false')
    parameters.put('selectedPolarisations', polarisation)
    parameters.put('outputSigmaBand', 'false')
    parameters.put('outputGammaBand', 'true')
    parameters.put('outputBetaBand', 'false')

    calibration = GPF.createProduct('Calibration', parameters, subset)
    print('Calibration done')

    ## Speckle-Filter

    #Variablen von oben hier benutzt.
    filterSizeX = '5'
    filterSizeY = '5'

    parameters = HashMap()

    parameters.put('sourceBands', 'Gamma0_{}'.format(polarisation))
    #parameters.put('filter', 'Lee')
    parameters.put('filter', 'Boxcar')
    parameters.put('filterSizeX', filterSizeX)
    parameters.put('filterSizeY', filterSizeY)
    parameters.put('dampingFactor', '2')
    parameters.put('estimateENL', 'true')
    parameters.put('enl', '1.0')
    parameters.put('numLooksStr', '1')
    parameters.put('targetWindowSizeStr', '3x3')
    parameters.put('sigmaStr', '0.9')
    parameters.put('anSize', '50')

    speckle_filter = GPF.createProduct('Speckle-Filter', parameters, calibration)
    print('Speckle-Filter done')

    ## Multilook
    azLooks = 3
    rgLooks = 3

    parameters = HashMap()
    parameters.put('sourceBands', 'Gamma0_{}'.format(polarisation))

    parameters.put('grSquarePixel', False)
    parameters.put('IndependentLooks', True)

    parameters.put('nRgLooks', rgLooks)
    parameters.put('nAzLooks', azLooks)
    parameters.put('outputIntensity', True)
    multilook = GPF.createProduct('Multilook', parameters, speckle_filter)
    print('Multilook done')

    ## Terrain-Correction
    parameters = HashMap()

    parameters.put('sourceBands', 'Gamma0_{}'.format(polarisation))
    parameters.put('demName', 'SRTM 3Sec')
    parameters.put('externalDEMFile', '')
    parameters.put('externalDEMNoDataValue', '0.0')
    parameters.put('externalDEMApplyEGM', 'true')
    parameters.put('demResamplingMethod', 'BILINEAR_INTERPOLATION')
    parameters.put('imgResamplingMethod', 'BILINEAR_INTERPOLATION')
    parameters.put('pixelSpacingInMeter', '30.0')
    parameters.put('pixelSpacingInDegree', '2.6949458523585647E-4')
    #parameters.put('pixelSpacingInDegree', '8.983152841195215E-5')
    parameters.put('mapProjection', 'AUTO:42001')
    parameters.put('nodataValueAtSea', 'true')
    parameters.put('saveDEM', 'false')
    parameters.put('saveLatLon', 'false')
    parameters.put('saveIncidenceAngleFromEllipsoid', 'false')
    parameters.put('saveProjectedLocalIncidenceAngle', 'false')
    parameters.put('saveSelectedSourceBand', 'true')
    parameters.put('outputComplex', 'false')
    parameters.put('applyRadiometricNormalization', 'false')
    parameters.put('saveSigmaNought', 'false')
    parameters.put('saveGammaNought', 'false')
    parameters.put('saveBetaNought', 'false')
    parameters.put('incidenceAngleForSigma0', 'Use projected local incidence angle from DEM')
    parameters.put('incidenceAngleForGamma0', 'Use projected local incidence angle from DEM')
    parameters.put('auxFile', 'Latest Auxiliary File')
    


    terrain_correction = GPF.createProduct('Terrain-Correction', parameters, multilook)
    print('Terrain-Correction done')

    ## Linear to dB
    parameters = HashMap()

    parameters.put('sourceBands', 'Gamma0_{}'.format(polarisation))
    parameters.put('outputImageScaleInDb', True)
    lintodB = GPF.createProduct('LinearToFromdB', parameters, terrain_correction)
    print('Linear to dB done')
    final_product = lintodB



    #return final_product
   
    ## Save the result

    #ProductIO.writeProduct(terrain_correction, s1_identifier + '.tif', 'GeoTIFF')
    #ProductIO.writeProduct(terrain_correction, s1_identifier + '.tif', "GeoTIFF-BigTIFF")
    #ProductIO.writeProduct(terrain_correction, s1_identifier + '.dim', 'BEAM-DIMAP')
    #print('Product saved as {}.dim'.format(s1_identifier))

    # ## Plot the result

    def plotBand(s1_identifier, product, rband, vmin, vmax):
        
        band = product.getBand(rband)

        w = band.getRasterWidth()
        h = band.getRasterHeight()

        band_data = np.zeros(w * h, np.float32)
        band.readPixels(0, 0, w, h, band_data)

        band_data.shape = h, w

        width = 12
        height = 12
        plt.figure(figsize=(width, height))
        imgplot = plt.imshow(band_data, cmap=plt.cm.binary, vmin=vmin, vmax=vmax)
        plt.axis('off')
        #plt.imshow(band_data, cmap=plt.cm.RdGy, vmin=vmin, vmax=vmax)
        
        plt.savefig(s1_identifier,bbox_inches='tight')
        #return imgplot 

    #plotBand(terrain_correction, 'Gamma0_{}'.format(polarisation), 0, 0.3)
    plotBand(s1_identifier, terrain_correction, 'Gamma0_{}'.format(polarisation), 0, 6)
    print('Plot saved as {}.png'.format(s1_identifier))
    return final_product #totally wrong hier move out plotter

def BandMath(stack):
    product = ProductIO.readProduct(stack)
    band_names = product.getBandNames()
    print(list(band_names))
    band1 = input("Enter Band 1:")
    band2 = input("Enter Band 2 which will be subtracted:")
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    ## Band Math
    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 2)
    #print('{}_minus_{}'.format(bandlist[x], bandlist[y]))
    targetBand1 = BandDescriptor()
    targetBand1.name = '{}_minus_{}'.format(band1, band2)
    targetBand1.type = 'float32'
    #targetBand1.expression = '(({} - {})<(-2))? 255 : 0'.format(band1, band2)
    targetBand1.expression = '(({} - {})<(-2))? 255 : 0'.format(band1, band2)
    parameters = HashMap()
    parameters.put('targetBands', targetBands)
    result = GPF.createProduct('BandMaths', parameters, product)
    print("Writing...")
    ProductIO.writeProduct(result, '{}_minus_{}.dim'.format(band1, band2), 'BEAM-DIMAP')
    print('{}_minus_{}.dim Done.'.format(band1, band2))

def BandMathList(stack):
    product = ProductIO.readProduct(stack)
    band_names = product.getBandNames()

    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

    BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

    bandlist = list(band_names)
    #targetBands = list()
    x=0
    y=1
    bandlength= len(bandlist)
    runs=(bandlength-1)
    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', runs)
    for i in bandlist:
        while y <= runs:
            #print('{}_minus_{}'.format(bandlist[x], bandlist[y]))
            #targetBand1 = '{}_minus_{}'.format(bandlist[x], bandlist[y])
            targetBand1 = BandDescriptor()
            targetBand1.name = '{}_minus_{}'.format(bandlist[x], bandlist[y])
            targetBand1.type = 'float32'
            targetBand1.expression = '(({} - {})<(-2))? 255 : 0'.format(bandlist[x], bandlist[y])
            
            print("Writing Band {} : {}_minus_{}".format(x,bandlist[x], bandlist[y]))
            #targetBands.append(targetBand1) 
            targetBands[x] = targetBand1
            
            x=x+1
            y=y+1

    """targetBand1 = BandDescriptor()
    targetBand1.name = 'first_{}_minus_last_{}'.format(bandlist[0], bandlist[bandlength])
    targetBand1.type = 'float32'
    targetBand1.expression = '(({} - {})<(-2))? 255 : 0'.format(bandlist[x], bandlist[y])
    print("Writing Band first_{}_minus_last_{}".format(bandlist[0], bandlist[bandlength])
    targetBands[bandlength] = targetBand1"""

    parameters = HashMap()
    parameters.put('targetBands', targetBands)
    result = GPF.createProduct('BandMaths', parameters, product)

    print("Writing...")
    ProductIO.writeProduct(result, 'BandMaths.dim', 'BEAM-DIMAP')
    print("BandMaths.dim Done.")
    ProductIO.writeProduct(result, 'BandMaths.tif', "GeoTIFF-BigTIFF")
    print("BandMaths.tif Done.")




def clean():
    # # Garbage collector
    thermal_noise_removal = None
    apply_orbit_file = None
    calibration = None
    speckle_filter = None
    terrain_correction = None
    product = None
    gc.collect()



if __name__ == "__main__":
    ziplist = zip_extract()
    #process(ziplist)
    #BandMathList(os.getcwd()+ "\create_stack.dim")
    BandMath(os.getcwd()+ "\create_stack.dim")
    
"""    for i in ziplist:
        i = os.path.splitext(i)[0]
        process(i)
    clean()"""