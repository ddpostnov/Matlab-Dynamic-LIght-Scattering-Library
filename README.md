Setting up the protocol
General notes:
•	Plan out the way to name files and folders in advance. It will save you time afterwards.  
Calibration recordings:
•	Framerate: 20Hz
•	Exposure time: adjust to reach 30% intensity saturation (<I>=90)
•	Duration: 5 min
•	Sample prep: piece of paper, preferably the same piece used and kept for the entire protocol.
•	Field of view: as large as possible, while ensuring homogeneous illumination
Imaging pulsatility
•	Framerate: 194Hz
•	Duration: 2 min or more
•	Animal prep: cranial window, awake is optimal
•	Field of view: 1024x512 (larger is possible when device throughput is disabled)
•	Notes: look out for dropped frames! Ensure physiological stability over the duration of the recording – otherwise you need to temporally crop the file.
Imaging vasomotion
•	Framerate: at least 25 frames per second, but the higher the better 
•	Duration: 10 min or more, 30 minutes is optimal
•	Animal prep: cranial window, awake is optimal
•	Field of view: any, but the larger the better
•	Note: having pulsatility recording prior to vasomotion is optimal to map arteries and veins
Imaging NVC
•	Framerate: at least 25 frames per second, but the higher the better
•	Duration: 120s baseline + 20x(5s 3Hz stim + 25 seconds recovery).
•	Animal prep: cranial window, awake is optimal but light ketamine is “okay”
•	Field of view: any, but the larger the better
•	Notes: Timing is important! Ensure some form of registration between stim timing and LSCI. If possible, longer recovery between stims and more repetitions are better. Make sure to monitor the animal condition, preferably record (and temporally register) pupil dilation and grooming behaviour.
Imaging slow dynamics in a single file (e.g. vasoreactivity)
•	Framerate: depends on the question, usually not less than 5
•	Duration: 30 minutes or more
•	Animal prep:  depends on question
•	Field of view: any, but the larger the better
•	Notes: if any stimuli (e.g. vasoreactivity) is used – ensure to register the timing with LSCI. Plan for other features (pulsatility vasomotion) and do additional baseline recordings as needed. Ensure position stability. If stability is impossible, it might be better to proceed with multiple files.

Imaging different conditions in multiple files (e.g. stroke)
•	Notes: Combination of the listed above. Ensure that the system configuration is maintained the same. Match field of view as close as possible. Prioritise tilt and focus matching over translation matching.
Longitudinal Imaging (multiple timepoints, e.g. ageing study) and absolute comparisons across between groups
•	Notes: Perform calibration recordings at chosen timepoints. Ideally on each experimental day. Ensure that the system configuration is maintained the same. Perform calibration recordings at chosen timepoints. Match field of view as close as possible. Prioritise tilt and focus matching over translation matching. 

 
Analysis
General notes:
•	All MATLAB files located in the folder “Dynamic Light Scattering Library v***” (DLSL) are required.
•	Make sure to delete all previous versions of the “Dynamic Light Scattering Library v***” before using a newer version. MATLAB might remember the path, and you will end up with conflicts. It is enough to have only one copy of the DLSL on PC as long as it accessible to all users.
•	As a user you will only need to interact with Launchers and Examples folders.
•	Make your own copy of required Launcher or Example file. You can store them anywhere you’d like. For your processing tasks – modify your Launcher and or Example files accordingly to the experimental protocol and your project needs. Leave the original Launchers and Examples unchanged.
•	Read comments (text after %) and, if relevant, headers (comments at the top) in scripts and functions. Do not use scripts blindly, get at least a basic understanding of what they are doing and what results you should expect.
•	In general – do not use multiple versions of the “Dynamic Light Scattering Library v***” when processing the files. Meaning if you want to use a new version – better to re-analyse the files from the first step. Compatibility between the versions is not guaranteed.
•	It helps to ask ChatGPT for advice if the code is difficult to comprehend. It is recommended to use ChatGPT o3 or other advanced versions that work well with coding. Note that if you use ChatGPT to generate code and it produces an error, sometimes you can just feed the error back to ChatGPT to produce a fix.
•	Addressing the files is often done using regular expressions (things like * and ? in the file name string). Read a bit about to understand how it is used. If you use it in combination with proper naming of the files – you can save a lot of time and reduce the chances using the wrong files. Again, ChatGPT is really good with regular expressions and can help you figuring the out.
•	When using processing pipelines, multiple .mat files will be generated. The  logic behind their names is that original file name is followed by a sequence of flags that are defined by the processing type and step. Examples:
o	_t_e_K_d.mat means (t) temporal contrast analysis; (e) estimated epoch; (K) contrast; (d) – file with 3D data.
o	_c_BFI_r.mat means (c) internal, usually cardiac, cycle; (BFI) blood flow index; (r) – file with processing results.
Generally, as a user, you would only interact with _BFI_r.mat files, while the processing steps would be defined by the project and the protocol.
•	At the final processing step, there is an option to export the key _BFI_r.mat results to an Excel table.
•	It is recommended to get “a feeling” of your data and some basic MATLAB skills. For that you can follow any processing pipeline to reach _BFI_d.mat file which will contain your 3D data (X,Y,Time) converted to contrast as well as the time vector. From there you can apply all sorts of data extraction and filtering procedures, such as ROI selection, or plot it as video or as an image. For inspiration check the Examples sourceOperations.mat or ask ChatGPT. An few exemplary prompts for ChatGPT are below:
o	“In MATLAB, I have a 3D matrix named source.data and a time vector named source.time. The source.data represents Blood Flow Index (BFI) obtained using laser speckle contrast imaging. The time is in seconds. I want to plot the image of that corresponds to the average of source.data in time. The colour limits should account for the possibility of the image containing some pixels with extremely high, extremely low or NaN values in respect to the meaningful data. Then I want to draw multiple regions of interest (ROI) interactively on the image. Finally I want all the ROI signals, which are averaged over pixels for individual regions, to be plotted versus time in two panels: first as a plot of actual values, second as a plot of normalised values, e.g. with each signal normalised from 0 to 1.”
o	“In MATLAB, I have a 3D matrix named source.data and a time vector named source.time. The source.data represents Blood Flow Index (BFI) obtained using laser speckle contrast imaging. Consider that some pixels might have abnormal values – very low (e.g. 0), very high (infinitiy or simply many times higher than all of the others) or NaNs. Make sure the further analysis corrects for it. The data consist of baseline period of 5 seconds, which is followed by a stimulation and response period. Find the peak response period and generate an image of response divided by the baseline. Highlight the area where response is above 10%.  Side by side make a plot of the average signal in the peak response area versus the remaining area, both normalised by respective baseline values. Explain how you found the peak of the response. “

Launchers
Launcher_pulsatility.m 
Pulsatility analysis provides steady state pulsatility and blood flow information. Currently it is designed only for data with rich vascular features and high quality images.
Are your recordings captured accordingly? 194 frames per second, cranial window etc?
If you record multiple conditions – are they all in individual files?
Ensure to copy the file to the location of choice and then edit it according to the project needs.

Launcher_slowDynamics.m. 
Either used to compare steady-state perfusion or responses to intervention.
Ensure to copy the file to the location of choice and then edit it according to the project needs.

Launcher_NVC.m. 
Used for analysis of externally triggered neurovascular coupling responses.
Ensure to copy the file to the location of choice and then edit it according to the project needs.

Launcher_basic.m 
Can be used for getting the feeling of your data without complex pre-processing. A minimal Launcher that consists of contrast calculation and BFI conversion. Won’t fit the requirements for most of the research projects.
Ensure to copy the file to the location of choice and then edit it according to the project needs.

Examples
Example_plotROIs.m 
Uses _d.mat file as an input (do not forget to set the file name) and _r.mat if available (for the categorical mask). Allows you to choose and plot ROIs. 

Example_GPT_plotROIs.m 
ChatGPT generated script for plotting ROIs. See the prompt in the previous section. 


Example_GPT_responseAnalysis.m 
ChatGPT generated script for analysing a response. See the prompt in the previous section. 

