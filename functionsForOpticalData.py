#======Imports========
import os
import numpy as np
import matplotlib.pyplot as plt

#======Functions========
#----Locating Experimental Data----
def locateData():
    """
    Function to auto-locate the exper data from the ~/InputData subDirectory
    Data expected in .txt format with 2 columns seperated by ";" in each row

    Returns
    -------
    numberOfSpectra : int
        Number of txt files found in subDir.
        
    spectrumDataList : list
        A list with the filePaths of all txt data found.
        
    txtFilesFoundSuccessfully: bool
        Returns True if .txt spectrum files were found and False if nothing was found
    """
   
    cwd = os.getcwd()  # locating current directory
    targDir = os.chdir(cwd+"/InputData") #change to subdirectory containing files
    
    #saving files to a list
    fileList = os.listdir(targDir) #List ALL the files in working directory
    
    #Filtering just .txt files from other types mistakenly being in this folder
    spectrumDataList=[]
    numberOfSpectra=0
    # filter for all the .txt files and store them
    for file in fileList:
        if file.endswith(".txt"):
            spectrumDataList.append(file)
            numberOfSpectra+=1 #spectra are counted
    
    #check if any .txt files were actually found
    txtFilesFoundSuccessfully=True
    if(len(spectrumDataList)==0):
        print("\nERROR:")
        print("No txt files with experimental data found!")
        print("Make sure your .txt spectra are located in subDir ~/InputData")
        txtFilesFoundSuccessfully=False
    return txtFilesFoundSuccessfully,numberOfSpectra,spectrumDataList
#end function


#---- Extracting data from first spectrum to act as a reference----
def TryReadFirstDataset(txtFileList):
    """
    Parameters
    ----------
    txtFileList : list
        List containing the filenames of the txt data found.

    Returns
    -------
    launchTime : list
        Contains the day,hour,minute and second of the first spectrum recorded.
        
    numberOfPointsInSpectrum : int
        The number of rows in each spectrum .txt file.
        
    Wavelegths : numpy.ndarray
        Contains values from 1st column in spectrum .txt files, 
        that corresponds to the wavelengths recorded by the optical sensor.
        
    IntensitiesRef : numpy.ndarray
        Contains values from 2nd column in the first "reference" spectrum .txt file, 
        that corresponds to the intensity recorded by the optical sensor.
    """
    
    #keep track of the rows in spectrum
    #This will be assumed to be a constant for every exper data file
    numberOfPointsInSpectrum=0
    
    #Waveleghts column is assumed the same in all files from the same experiment
    Wavelegths,IntensitiesRef=[],[]
    try:
        #record launch time from first spectrum file name
        launchTime = txtFileList[0].split("_")
        launchTime[0] = launchTime[0][6:]
        #--------------------
        file = open(txtFileList[0], 'rt')
        for myline in file:
            numberOfPointsInSpectrum+=1
            #txt has a specific format with 2 columns seperated by ;
            wavelegth,intensity=myline.rstrip('\n').split(";")
            #extract the numbers, format them properly and store them
            Wavelegths.append( float(wavelegth))
            IntensitiesRef.append( float(intensity))
    finally:
        file.close()
    
    #change to numpy arrays
    Wavelegths = np.array(Wavelegths)
    IntensitiesRef=np.array(IntensitiesRef)
    return launchTime,numberOfPointsInSpectrum,Wavelegths,IntensitiesRef
#end function

def wavelengthToIndex(reqWav,Wavelengths):
    """
    To be used by the functions module
    DO NOT CALL
    """
    searchIndex=0
    while(Wavelengths[searchIndex]<reqWav): searchIndex+=1
    return searchIndex
#end function

#Visualize the first Spectrum
def PlotFirstSpectrum(Wavelengths,IntensitiesRef):
    """
    Plots the Intestity vs Wavelength in the first spectrum recorded
    
    Parameters
    ----------
    Wavelegths : numpy.ndarray
        Contains values from 1st column in spectrum .txt files, 
        that corresponds to the wavelengths recorded by the optical sensor.
        
    IntensitiesRef : numpy.ndarray
        Contains values from 2nd column in the first "reference" spectrum .txt file, 
        that corresponds to the intensity recorded by the optical sensor.

    Returns
    -------
    None.

    """
    fig,ax=plt.subplots()
    ax.set_title("First Optical Spectrum")
    ax.set_xlabel("Wavelength (nm)");ax.set_ylabel('Intensity')
    plt.plot(Wavelengths,IntensitiesRef,ms=0.5)
#end function

 
#Function assisting time-related calculations (Not to be called by the user)
def calculateTimeSinceLaunch(launchTime,time):
    """
    To be used by the functions module
    DO NOT CALL
    """
    return 24*(float(time[0])-float(launchTime[0]))+ (float(time[1])-float(launchTime[1]))+(float(time[2])-float(launchTime[2]))/60 +(float(time[3])-float(launchTime[3]))/3600
#end function
    

#----Extract all intensity data----
def ScanFiles(txtFileList,numberOfSpectra,numberOfPointsInSpectrum,launchTime,Wavelengths,integrStartWav,integrFinishWav,selectedWavelengths):
    """
    Function to extract data from all spectra
        
    Parameters
    ----------
    txtFileList : list
        A list with the filePaths of all txt data found.
        
    numberOfSpectra : int
        Number of txt files found in subDir.
    
    numberOfPointsInSpectrum : int
        The number of rows in each spectrum .txt file. 

    launchTime : list
        Contains the day,hour,minute and second of the first spectrum recorded.
        
    Wavelengths : numpy.ndarray
        Contains values from 1st column in spectrum .txt files, 
        that corresponds to the wavelengths recorded by the optical sensor.
        
    integrStartWav : float
        First wavelength in ROI for averaging.
        
    integrFinishWav : float
        Last wavelength in ROI for averaging.
        
    selectedWavelengths : list
        Contains wavelengths of interest.

    Returns
    -------
    timesFromLaunchInHours : numpy.ndarray
        Contains the time (in hours) passed since launch for each spectrum.
        
    intensities : numpy.ndarray (2D)
        Contains values from 2nd column in the EVERY spectrum .txt file, 
        that corresponds to the intensity recorded by the optical sensor.
        Each Row is one spectrum
        
    averageIntensities : numpy.ndarray
        Contains the intesity averaged in ROI for each spectrum .
    """
    
    #-allocate memory to store the data
    timesFromLaunchInHours=np.zeros(numberOfSpectra)
    intensities = np.zeros( (numberOfSpectra,numberOfPointsInSpectrum),dtype=np.float64)
    averageIntensities    =np.zeros(numberOfSpectra)
    
    #-find useful indices
    specIndex=0
    firstIndex=wavelengthToIndex(integrStartWav,Wavelengths)
    lastIndex =wavelengthToIndex(integrFinishWav,Wavelengths)
    selectedIndices=[]
    for selectedWav in selectedWavelengths:
        selectedIndices.append(wavelengthToIndex(selectedWav,Wavelengths))        

    integrFinishWav=1000
    for txtFile in txtFileList: #For Each Spectrum
        #get time from filenames    
        time = txtFile.split("_")
        time[0] = time[0][6:]
        #time calculation 
        timesFromLaunchInHours[specIndex]=calculateTimeSinceLaunch(launchTime,time)
        #---------------------------------------
        arrIndex=0
        try:
            file = open(txtFile, 'rt')
            for line in file:  
                #extract the data
                wavelength,intensity=line.rstrip('\n').split(";")
                intensities[specIndex,arrIndex]= intensity
                arrIndex+=1
            #end for (lines)
            
        finally:
            file.close()
            
        #get the average intensity 
        averageIntensities[specIndex]= np.mean(intensities[specIndex,firstIndex:lastIndex])
        specIndex+=1
        
        if(specIndex%50 == 0):print(f"{specIndex}/{len(txtFileList)} files scanned")
    #end for (files)
    print(f"all txt files scanned",end="\n\n")
    return  timesFromLaunchInHours,intensities,averageIntensities            
#end function

#-------------------------------------
def plotAverageIntensityOverTime(timesFromLaunchInHours,averageIntensities,integrStartWav,integrFinishWav):
    """
    Parameters
    ----------
    timesFromLaunchInHours : numpy.ndarray
        Contains the time (in hours) passed since launch for each spectrum.
    averageIntensities : numpy.ndarray
        Contains the intesity averaged in ROI for each spectrum .

    Returns
    -------
    None.

    """
    fig2,ax2=plt.subplots()
    ax2.set_xlabel("time (h)");ax2.set_ylabel('Average Intensity')
    ax2.scatter(timesFromLaunchInHours,averageIntensities,s=1,label=f"ROI:{integrStartWav}-{integrFinishWav} nm")
    ax2.legend()
    plt.show()    
#end function

#----------------------------------------
def getSelectedWavelengthsAtTimeStamp(requestedElapsedTime,timesFromLaunchInHours,selectedWavelengths,wavelengths,intensities):
    """
    To be used internally in the functions module
    DO NOT CALL
    """
    #locating (by index) the spectrum with time closest to the requested time
    timeIndex=0
    while(timesFromLaunchInHours[timeIndex]<requestedElapsedTime): 
        timeIndex+=1
    #return the intensities of the selected Wavelengths for this spectrum
    selectedIntensities=[]
    for selectedWavelength in selectedWavelengths:
        selectedIntensity = intensities[timeIndex,wavelengthToIndex(selectedWavelength,wavelengths)]
        selectedIntensities.append(selectedIntensity)
    #endFor
    return selectedIntensities
#end function

def printUnloadedVsLoadedComparison(unloadedTimeStamp,loadedTimeStamp,timesFromLaunchInHours,selectedWavelengths,wavelengths,intensities):
    """
    DEPRECATED
    DO NOT USE
    """
    unloadedSelectedIntensities = getSelectedWavelengthsAtTimeStamp(unloadedTimeStamp,timesFromLaunchInHours,selectedWavelengths,wavelengths,intensities)
    loadedSelectedIntensities   = getSelectedWavelengthsAtTimeStamp(  loadedTimeStamp,timesFromLaunchInHours,selectedWavelengths,wavelengths,intensities)
    print()
    print(f"Unloaded VS Loaded Sample for Selected Wavelengths")
    for i in range(len(unloadedSelectedIntensities)):
        ratio = loadedSelectedIntensities[i]/unloadedSelectedIntensities[i]
        lnRatio= -np.log(loadedSelectedIntensities[i]/unloadedSelectedIntensities[i])
        print(f"wv: {selectedWavelengths[i]} nm: L/UnL: {ratio  } -ln(I(c)/Io) = { lnRatio }")
#end function

#-----------------------------------------------------------------------------
def GetIntensitiesAtTimestamp(requestedElapsedTime,timesFromLaunchInHours,intensities,spectrumDataList):
    """
    To be used internally in the functions module
    DO NOT CALL
    """
    #locating (by index) the spectrum with time closest to the requested time
    timeIndex=0
    while(timesFromLaunchInHours[timeIndex]<requestedElapsedTime): 
        timeIndex+=1
    print(f"Requested timeStamp {requestedElapsedTime} h,  file found: {spectrumDataList[timeIndex]} ")
    return intensities[timeIndex,:]
#end function-----------------------------------------------------------------


def GetStatisticsNearSelectedWavelengths(wavelengths,intensitiesSpec,selectedWavelengths):
    """
    To be used internally in the functions module
    DO NOT CALL
    """
    offset=10
    meanIntensities=[]
    stdsIntensities=[]
    print("Intensity statistics near the selected wavelegths for this timestamp:")
    for selectedWavelength in selectedWavelengths:
        #find the indices corresponding to wavls near the selected ones
        startIndex  = wavelengthToIndex(selectedWavelength, wavelengths) - offset
        finishIndex = wavelengthToIndex(selectedWavelength, wavelengths) + offset
        #
        intensitiesNearSelected= np.array( intensitiesSpec[startIndex:finishIndex])
        m = np.mean(intensitiesNearSelected)
        s = np.std(intensitiesNearSelected)
        print(f"wv: {selectedWavelength} nm mean: {m} std: {s}")
        meanIntensities.append(m)
        stdsIntensities.append(s)
    print()
    return meanIntensities, stdsIntensities


def ComapreStatisticsOfSelectedWavelengths(spectrumDataList,timesFromLaunchInHours,wavelengths,intensities,selectedWavelengths,unloadedTimeStamp,loadedTimeStamp):
    """

    Parameters
    ----------
    spectrumDataList : list
        A list with the filePaths of all txt data found.
        
    timesFromLaunchInHours : numpy.ndarray
        Contains the time (in hours) passed since launch for each spectrum.
    
    wavelengths : 
        Contains values from 1st column in spectrum .txt files, 
        that corresponds to the wavelengths recorded by the optical sensor.
    
    intensities : numpy.ndarray (2D)
        Contains values from 2nd column in the EVERY spectrum .txt file, 
        that corresponds to the intensity recorded by the optical sensor.
        Each Row is one spectrum
    
    selectedWavelengths : list
        Contains wavelengths of interest.
        
    unloadedTimeStamp : float
        TimeStamp corresponding to the sample before it is loaded with H2.
    
    loadedTimeStamp : float
        TimeStamp corresponding to the sample before it is loaded with H2.
    

    Returns
    -------
    None.

    """
    #get statistics from the "unloaded" timestamp 
    temparr1 = GetIntensitiesAtTimestamp(unloadedTimeStamp,timesFromLaunchInHours,intensities,spectrumDataList)
    m_un,s_un = GetStatisticsNearSelectedWavelengths(wavelengths,temparr1,selectedWavelengths)
    
    #get statistics from the "loaded" timestamp 
    temparr2 = GetIntensitiesAtTimestamp(loadedTimeStamp  ,timesFromLaunchInHours,intensities,spectrumDataList)
    m_l,s_l =GetStatisticsNearSelectedWavelengths(wavelengths,temparr2,selectedWavelengths)
    
    print(f"Unloaded VS Loaded Sample for Selected Wavelengths")
    for i in range(len(m_un)):
        ratio = m_l[i]/m_un[i]
        lnRatio= -np.log(m_l[i]/m_un[i])
        
        #error propagation
        err_ratio = ratio * np.sqrt( (s_l[i]/m_l[i])**2 +(s_un[i]/m_un[i])**2  )
        err_lnRatio = err_ratio/ratio
        print(f"wv: {selectedWavelengths[i]} nm: L/UnL: {ratio }, err {err_ratio}  -ln(I(c)/Io) = { lnRatio } err:{err_lnRatio}")
    
    return

