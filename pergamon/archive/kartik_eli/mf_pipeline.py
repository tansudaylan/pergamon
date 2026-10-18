from tdpy.verbosity import print
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math
import mock_light_curves as mlc
import confusion_matrix as cf
import os
import scipy.signal
import time
import matplotlib.backends.backend_pdf as pdfs

class Light_Curve():
    kernel = 99
    alpha = 0.1
    beta = 0.25
    centroid_threshold = 0.5
    
    def __init__(self, tess_bjds, pdcsap_fluxes, centroid, templates):
        self.tess_bjds = tess_bjds
        self.pdcsap_fluxes = pdcsap_fluxes
        self.centroid = centroid
        
        self.flat_lcur = []
        self.valid_times = []
        self.valid_fluxes = []
        self.valid_centroid = []
        self.flat_centroid = []
        
        self.start_edge = 0
        self.end_edge = 0
        self.left_edge = 0
        self.right_edge = 0
        
        self.templates = templates
        self.result = False
        self.correlations = []
        self.best_correlation = 0
        self.best_template = self.templates[0]
        self.window = 0
        
        self.flag = "Normal"
        self.detection_indices = []
        self.real_detections = []
        self.ambiguous_real = []
        self.flare_detections = []
        self.ambiguous_flare = []
        self.centroid_detections = []
        
    
    def get_valid_data(self):
        notnan_times = []
        notnan_fluxes = []
        notnan_centroid = []
        
        for i in range(len(self.pdcsap_fluxes)):
            flux = self.pdcsap_fluxes[i]
            if not np.isnan(flux):
                notnan_times.append(self.tess_bjds[i])
                notnan_fluxes.append(self.pdcsap_fluxes[i])
                notnan_centroid.append(self.centroid[i])
        
        return notnan_times, notnan_fluxes, notnan_centroid
                            
    def sigma_clip(self, flat_lcur_unclipped, notnan_times, notnan_fluxes, notnan_centroid):
        sigma_clip = 10*np.std(flat_lcur_unclipped)
        length = len(flat_lcur_unclipped)
        flat_lcur = []
        valid_times = []
        valid_fluxes = []
        valid_centroid = []
        for i in range(length):
            flux = flat_lcur_unclipped[i]
            if abs(flux) > sigma_clip:
                if i == 0:
                    if flat_lcur_unclipped[i+1] > sigma_clip:
                        flat_lcur.append(flux)
                        valid_times.append(notnan_times[i])
                        valid_fluxes.append(notnan_fluxes[i])
                        valid_centroid.append(notnan_centroid[i])
                elif i == length-1:
                    if flat_lcur_unclipped[i-1] > sigma_clip:
                        flat_lcur.append(flux)
                        valid_times.append(notnan_times[i])
                        valid_fluxes.append(notnan_fluxes[i])
                        valid_centroid.append(notnan_centroid[i])
                else:
                    if flat_lcur_unclipped[i+1] > sigma_clip or flat_lcur_unclipped[i-1] > sigma_clip:
                        flat_lcur.append(flux)
                        valid_times.append(notnan_times[i])
                        valid_fluxes.append(notnan_fluxes[i])
                        valid_centroid.append(notnan_centroid[i])
            else:
                flat_lcur.append(flux)
                valid_times.append(notnan_times[i])
                valid_fluxes.append(notnan_fluxes[i])
                valid_centroid.append(notnan_centroid[i])
        
        return flat_lcur, valid_times, valid_fluxes, valid_centroid

    def get_edges(self, notnan_times):
        left_edge = None
        right_edge = None
        for k in range(len(notnan_times)):
            if notnan_times[k] - notnan_times[k-1] > 1:
                left_edge = notnan_times[k-1]
                right_edge = notnan_times[k]
                return left_edge, right_edge
            
    def match_filter(self):
        initial = True
        prev_positive_corr = None
        prev_result = False
        
        # perform cross-correlation for all template widths
        for template in self.templates:  
            correlations = scipy.signal.correlate(self.flat_lcur, template, mode='valid')
            highest_corr = max(correlations)
            result = highest_corr > Light_Curve.alpha
        
            # choose best correlation result so far 
            if initial or highest_corr > self.best_correlation:
                self.best_correlation = highest_corr
                self.result = result
                self.correlations = correlations
                self.best_template = template
                initial = False
                
            #break if new template has lower correlation than previous template
            if prev_result:
                if highest_corr < prev_positive_corr:
                    break
            if result:
                prev_positive_corr = highest_corr
                
            prev_result = result
            
    def get_detection_indices(self):
        detection_indices = []
        prev_detection = False
        best_detection_info = None
        need_to_add = False
        corr_length = len(self.correlations)
        
        #find locations of positive detections
        for j in range(corr_length):
            if self.correlations[j] > Light_Curve.alpha:
                if not prev_detection:
                    need_to_add = True
                    best_detection_info = (j, self.correlations[j])
                else:
                    if self.correlations[j] > best_detection_info[1]:
                        best_detection_info = (j, self.correlations[j])
                prev_detection = True
            else:
                if need_to_add:
                    detection_indices.append(best_detection_info[0])
                    need_to_add = False
                prev_detection = False
        
        return detection_indices
    
    def classify_detections(self):
        nonedge_detections = []
        # find locations of non-edge detections
        for detection in self.detection_indices:
            t = self.valid_times[detection]
            if not self.is_edge_detection(t, self.start_edge, self.end_edge, self.left_edge, self.right_edge):
                nonedge_detections.append(detection)
                
        # all detections are edge detections so we flag as 'edge'
        if len(nonedge_detections) == 0:
            self.flag = 'edge'
            
        min_corr = min(self.correlations)
        # flag as transit if there is a stronger negative correlation than a positive one
        if self.flag == 'SL':
            if abs(min_corr) > 1.5 * self.best_correlation:
                self.flag = 'transit'
                
        # check centroid of detection against template to see if this is a valid detection
        candidate_detections = []
        if self.flag == 'SL':
            for detection in nonedge_detections:
                centroid_sample = mlc.offset_and_normalize(self.flat_centroid[detection:detection+self.window])
                centroid_template = mlc.offset_and_normalize(mlc.generate_centroid_template(1, self.window//3))
                centroid_template2 = np.flip(centroid_template)
                #detection_sample = mlc.offset_and_normalize(flat_lcur[detection:detection+window])
                centroid_corr = scipy.signal.correlate(centroid_sample, centroid_template, mode='valid')[0]
                centroid_corr2 = scipy.signal.correlate(centroid_sample, centroid_template2, mode='valid')[0]
                if centroid_corr >= Light_Curve.centroid_threshold or centroid_corr2 >= Light_Curve.centroid_threshold:
                    self.centroid_detections.append((detection, self.best_template, None, None))
                else:
                    candidate_detections.append(detection)
            # flag as centroid if detections are centroid detections
            if len(candidate_detections) == 0:
                self.flag = 'centroid'
                
        # find locations of all real and flare detections
        ambiguity = .05
                    
        for detection in candidate_detections:
            time_window = self.valid_times[detection:detection+self.window]
            detection_sample = mlc.offset_and_normalize(self.flat_lcur[detection:detection+self.window])
            gaussian_template = mlc.offset_and_normalize(mlc.generate_template(1, self.window//3, gaps=True, detection_time=self.valid_times[detection+self.window//2], times=time_window))
            #flare_template = mlc.generate_flare_template(1, window//3, gaps=True, detection_time=valid_times[detection+window//2], times=time_window)
            gaussian_result = scipy.signal.correlate(detection_sample, self.best_template, mode='valid')[0]
            flare_result = 0
            best_flare_template = None
            best_flare_time_window = None
            for i in range(-5, 6):
                flare_time_window = self.valid_times[detection+i:detection+self.window+i]
                flare_template = mlc.offset_and_normalize(mlc.generate_flare_template(1, self.window//3, gaps=True, detection_time=self.valid_times[detection+self.window//2+i], times=flare_time_window))
                flare_window = mlc.offset_and_normalize(self.flat_lcur[detection+i:detection+self.window+i])
                flare_correlation = scipy.signal.correlate(flare_window, flare_template, mode='valid')[0]
                if flare_correlation > flare_result:
                    flare_result = flare_correlation
                    best_flare_template = flare_template
                    best_flare_time_window = flare_time_window
            # compare correlation to gaussian and flare templates
            if flare_result > gaussian_result:
                # check for ambiguity between the gaussian and flare results
                if abs(flare_result - gaussian_result) < ambiguity:
                    self.ambiguous_flare.append((detection, gaussian_template, best_flare_template, best_flare_time_window))
                else:
                    self.flare_detections.append((detection, gaussian_template, best_flare_template, best_flare_time_window))
            else:
                if abs(gaussian_result - flare_result) < ambiguity:
                    self.ambiguous_real.append((detection, gaussian_template, best_flare_template, best_flare_time_window))
                else:
                    self.real_detections.append((detection, gaussian_template, best_flare_template, best_flare_time_window))
                            
        # all detections are not SL detections so we flag as 'flare', 'ambiguousFlare', 'ambiguousSL'                
        if self.flag == 'SL' and len(self.real_detections) == 0:
            if len(self.ambiguous_real) > 0:
                self.flag = 'ambiguousSL'
            elif len(self.ambiguous_flare) > 0:
                self.flag = 'ambiguousFlare'
            else:
                self.flag = 'flare'
                
        # check if the detections exceed the higher threshold
        if self.flag == 'SL':
            for detection in self.real_detections:
                detection_corr = self.correlations[detection[0]]
                if detection_corr > Light_Curve.beta:
                    self.flag = 'highSL'
                    break
    
    def is_edge_detection(self, t, start, end, left, right, cutoff=.5):
        return any([t <= start + cutoff, 
                    t >= end - cutoff, 
                    (t >= left - cutoff) and (t <= left),
                    (t <= right + cutoff) and (t >= right)])
        
    def run_pipeline(self):
        # Take only valid data points       
        notnan_times, notnan_fluxes, notnan_centroid = self.get_valid_data()
                
        # determine left and right edges for orbital gap
        self.left_edge, self.right_edge = self.get_edges(notnan_times)
                
        # subtract median filter to get high frequency noise for light curve
        flat_lcur_unclipped = mlc.offset_and_normalize(notnan_fluxes - scipy.signal.medfilt(notnan_fluxes, Light_Curve.kernel))
        
        # remove lone outliers to avoid fitting to outliers
        self.flat_lcur, self.valid_times, self.valid_fluxes, self.valid_centroid = self.sigma_clip(flat_lcur_unclipped, notnan_times, notnan_fluxes, notnan_centroid)
        
        # subtract median filder to get high frequency noise for centroid
        self.flat_centroid = mlc.offset_and_normalize(self.valid_centroid - scipy.signal.medfilt(self.valid_centroid, Light_Curve.kernel))
          
        # run match filter on the flat light curve with all templates
        self.match_filter()
        
        if self.result:
            self.window = len(self.best_template)
            self.start_edge = self.valid_times[0]
            self.end_edge = self.valid_times[len(self.correlations)-1]
            self.flag = 'SL'
            
            self.detection_indices = self.get_detection_indices()
            self.classify_detections()
            
            
        
def mf_pipeline(directory, result_foldername, mock=False, num_simulations=None):
    '''
    Pipeline runs a match filter on light curves and finds if the light curve 
    matches a predetermined template.
    
    directory: string, directory where light curve .fits files are located
    result_foldername: string, location of resulting plots, data 
    mock: bool, true if using mock data, False if using real light curve
    num_simulations: int, number of mock light curves to generate (if mock)
        
    returns: dict, results for each light curve (file)
    '''
    kernel = 99
    # generate templates varying width from 30mins to 2hrs (15 bins to 60 bins)
    templates = []
    widths = [5*j + 15 for j in range(10)] + [10*j + 70 for j in range(4)]
    for width in widths: 
        templates.append(mlc.offset_and_normalize(mlc.generate_template(1, width)))
        
    if mock:
        num_bins = 10
        counter = 0
        
        # create folder for resulting plots
        if not os.path.exists(result_foldername):
            os.mkdir(result_foldername)
        
        # results array for completeness analysis of each variable
        i_actual = [[] for _ in range(num_bins)]
        i_predicted = [[] for _ in range(num_bins)]
        inclinations = [[] for _ in range(num_bins)]
        
        P_actual = [[] for _ in range(num_bins)]
        P_predicted = [[] for _ in range(num_bins)]
        periods = [[] for _ in range(num_bins)]
        
        mbh_actual = [[] for _ in range(num_bins)]
        mbh_predicted = [[] for _ in range(num_bins)]
        mbhs = [[] for _ in range(num_bins)]
        
        ms_actual = [[] for _ in range(num_bins)]
        ms_predicted = [[] for _ in range(num_bins)]
        mss = [[] for _ in range(num_bins)]
        
        noise_actual = [[] for _ in range(num_bins)]
        noise_predicted = [[] for _ in range(num_bins)]
        noises = [[] for _ in range(num_bins)]
        
        threshold_actual = [[] for _ in range(num_bins)]
        threshold_predicted = [[] for _ in range(num_bins)]
        thresholds = [[] for _ in range(num_bins)]
        
        total_actual = []
        total_predicted = []
        
        # generate bins
        cosi_bins = [j/num_bins * .01  for j in range(num_bins)]
        P_bins = [np.e**(j*(np.log(27)-np.log(1))/num_bins) for j in range(num_bins)]
        mbh_bins = [15*j/num_bins + 5 for j in range(num_bins)]
        ms_bins = [j/num_bins + 0.5 for j in range(num_bins)]
        noise_bins = np.array([0,50,100,150,250,500,750,1000,1500,2000]) 
        altered_noise_bins = noise_bins * 10**(-6) + 10**(-9)
        threshold_bins = np.array([.01*i for i in range(10)]) 
        
        start = time.time()
        for z in range(1, num_simulations+1):
            
            # randomly generate relevant parameters from known priors
            P = mlc.P_rng()
            M_BH = mlc.mbh_rng()
            M_S = mlc.ms_rng()
            i = mlc.i_rng()
            cosi = math.cos(i)
            noise = np.random.choice(altered_noise_bins)
            threshold = np.random.random() * .1
            
            # generate positive/negative signal at random
            pos_signal = np.random.choice([True, False])
            if pos_signal:
                lcur, EV, Beam, SL = mlc.generate_light_curve(P, i, M_BH, M_S=M_S, std=noise)
            else:
                lcur = mlc.generate_flat_signal(noise)
                
            # subtract median filter from signal and normalize for correlation analysis
            flat_lcur = lcur - scipy.signal.medfilt(lcur, kernel)
            flat_lcur = mlc.offset_and_normalize(flat_lcur)
            
            initial = True
            best_result = None
            best_corr = 0
            best_correlations = None
            best_template = None
            
            # perform cross-correlation for all template widths
            for template in templates:  
                correlations = scipy.signal.correlate(flat_lcur, template)
                highest_corr = max(correlations)
                result  = highest_corr > threshold
                    
                # choose best correlation result so far (break on positive signal detection)
                if initial or highest_corr > best_corr:
                    best_corr = highest_corr
                    best_result = result
                    best_correlations = correlations
                    best_template = template
                    initial = False
                    
                if result:
                    break
            
            # plot some light curves and their correlations at random based on prediction
            if (best_result != pos_signal and np.random.random() > .98) or (pos_signal and np.random.random() > .99):
                lc_folder = "./{}/lc{}".format(result_foldername, counter)
                if not os.path.exists(lc_folder):
                    os.mkdir(lc_folder)
                mlc.plot_lc(lcur, P, M_BH, i, M_S, filename="{}/lcur{}.pdf".format(lc_folder, counter), EV=EV if pos_signal else None, Beam=Beam if pos_signal else None, SL=SL if pos_signal else None)
                mlc.plot_lc(flat_lcur, P, M_BH, i, M_S, filename="{}/flat_lcur{}.pdf".format(lc_folder, counter))
                mlc.plot_corr(best_correlations, P, M_BH, i, M_S, threshold, noise, "{}/corr{}.pdf".format(lc_folder, counter))
                counter += 1
            
            # bin result depending on parameter values
            if pos_signal or best_result:
                i_binned, P_binned, mbh_binned, ms_binned = False, False, False, False
                for k in range(num_bins-1,-1,-1):
                    if not i_binned and cosi >= cosi_bins[k]:
                        inclinations[k].append(i)
                        i_actual[k].append(pos_signal)
                        i_predicted[k].append(result)
                    if not P_binned and P >= P_bins[k]:
                        periods[k].append(i)
                        P_actual[k].append(pos_signal)
                        P_predicted[k].append(result)
                    if not mbh_binned and M_BH >= mbh_bins[k]:
                        mbhs[k].append(i)
                        mbh_actual[k].append(pos_signal)
                        mbh_predicted[k].append(result)
                    if not ms_binned and M_S >= ms_bins[k]:
                        mss[k].append(i)
                        ms_actual[k].append(pos_signal)
                        ms_predicted[k].append(result)
                    if all([i_binned, P_binned, mbh_binned, ms_binned]):
                        break
                    
            threshold_binned, noise_binned = False, False
            for k in range(num_bins-1,-1,-1):
                if not threshold_binned and threshold >= threshold_bins[k]:
                    thresholds[k].append(i)
                    threshold_actual[k].append(pos_signal)
                    threshold_predicted[k].append(result)
                if not noise_binned and noise == altered_noise_bins[k]:
                    noises[k].append(i)
                    noise_actual[k].append(pos_signal)
                    noise_predicted[k].append(result)
                if noise_binned and threshold_binned:
                    break
                    
            total_actual.append(pos_signal)
            total_predicted.append(best_result) 
            
            # perform completeness analysis
            if z%1000 == 0 and z != 0:
                print('{} simulations complete'.format(z))
                prefix = "./{}/".format(result_foldername)
                plot_completeness(r'$\alpha$', threshold_bins, threshold_actual, threshold_predicted, num_bins, prefix + 'alpha', scale=True)
                plot_completeness('Noise [ppm]', noise_bins, noise_actual, noise_predicted, num_bins, prefix + 'noise', scale=True)
                plot_completeness('cosi', cosi_bins, i_actual, i_predicted, num_bins, prefix + 'cosi', scale=True)
                plot_completeness('Period [days]', P_bins, P_actual, P_predicted, num_bins, prefix + 'period')
                plot_completeness(r'$M_{BH} [M_{\odot}]$', mbh_bins, mbh_actual, mbh_predicted, num_bins, prefix + 'mbh')
                plot_completeness(r'$M_{\star} [M_{\odot}]$', ms_bins, ms_actual, ms_predicted, num_bins, prefix + 'ms')
            
        end = time.time()
        print("{} minutes".format(round((end - start)/60, 2)))
        
        return total_actual, total_predicted
            
    else:
        counter = 1
        num_files = 0
        flags = ['SL', 'edge', 'transit', 'flare', 'highSL', 'ambiguousSL', 'ambiguousFlare', 'centroid']
        results = {flag: set() for flag in flags}
        if not os.path.exists(directory+result_foldername):
            os.mkdir(directory+result_foldername)
            for flag in flags:
                os.mkdir(directory + result_foldername + '/' + flag)
            
        for filename in os.listdir(directory):
            num_files += 1
            if num_files%500 == 0:
                print('{} files completed'.format(num_files))
                
            if filename.endswith(".fits"):
                fits_file = directory + filename
                try:
                    with fits.open(fits_file, mode="readonly", memmap=False) as hdulist:
                        tess_bjds = hdulist[1].data['TIME']
                        #sap_fluxes = hdulist[1].data['SAP_FLUX']
                        pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
                        centroid = hdulist[1].data['MOM_CENTR1']
                except:
                    print('Could not open file: {}'.format(filename))
                    continue
                
                # Create a light curve object with the file data and run pipeline on it
                lcur_object = Light_Curve(tess_bjds, pdcsap_fluxes, centroid, templates)
                lcur_object.run_pipeline()
                            
                if lcur_object.result:
                    results[flag].add(filename)
                    folder = '{}{}/{}'.format(directory, result_foldername, flag)
                    pdf = pdfs.PdfPages('{}/light_curve{}.pdf'.format(folder, counter))
                    counter += 1
                         
                    # create plot for full light curve data
                    fig, ax = plt.subplots(5, sharex=True, figsize=(6, 10))
                    #fig.suptitle(filename)
                    
                    # plot light curve
                    ax[0].set_title(filename)
                    ax[0].plot(lcur_object.valid_times, lcur_object.valid_fluxes, 'ko', rasterized=True, markersize=1)
                    ax[0].set_ylabel('PDCSAP Flux')
                    
                    # plot flat light curve
                    ax[1].plot(lcur_object.valid_times, lcur_object.flat_lcur, 'ko', rasterized=True, markersize=1)
                    ax[1].set_ylabel('Relative Flux')
                    
                    # plot correlation
                    ax[2].plot(lcur_object.valid_times[:len(lcur_object.correlations)], lcur_object.correlations, 'ko', rasterized=True, markersize=1)
                    ax[2].plot([lcur_object.valid_times[0], lcur_object.valid_times[len(lcur_object.correlations)-1]], [Light_Curve.alpha, Light_Curve.alpha], '--', color='orange', rasterized=True)
                    if lcur_object.flag == 'highSL':
                        ax[2].plot([lcur_object.valid_times[0], lcur_object.valid_times[len(lcur_object.correlations)-1]], [Light_Curve.beta, Light_Curve.beta], 'b--', rasterized=True)
                    ax[2].set_ylabel('Correlation')
                      
                    # plot centroid
                    ax[3].plot(lcur_object.valid_times, lcur_object.valid_centroid, 'ko', rasterized=True, markersize=1)
                    ax[3].set_ylabel('Centroid')
                    
                    # plot flat centroid
                    ax[4].plot(lcur_object.valid_times, lcur_object.flat_centroid, 'ko', rasterized=True, markersize=1)
                    ax[4].set_ylabel('Relative Centroid')
                    ax[4].set_xlabel('Time [days]')
                    
                    plt.tight_layout()
                    pdf.savefig(fig)
                    plt.close()
                    
                    # zoomed in plot on location of each real positive detection
                    for detection in lcur_object.real_detections:
                        plot_detection(detection[0], lcur_object.window, 'SL', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each ambiguous positive detection
                    for detection in lcur_object.ambiguous_real:
                        plot_detection(detection[0], lcur_object.window, 'Ambiguous SL', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each ambiguous flare detection
                    for detection in lcur_object.ambiguous_flare:
                        plot_detection(detection[0], lcur_object.window, 'Ambiguous Flare', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each flare detection
                    for detection in lcur_object.flare_detections:
                        plot_detection(detection[0], lcur_object.window, 'Flare', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                    # zoomed in plot on location of each centroid detection
                    for detection in lcur_object.centroid_detections:
                        plot_detection(detection[0], lcur_object.window, 'Centroid', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.best_correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                    
                    pdf.close()
                    
        # make a pie chart of the distribution of light curves in each bin
        pie_slices = [len(results[flag]) for flag in flags]
        plt.figure()
        plt.title('Distribution of Positive Detections')
        plt.pie(pie_slices, labels=flags)
        plt.savefig(directory + result_foldername + "/distribution.pdf")
        plt.close()
    

def plot_completeness(variable, bins, actual, predicted, num_bins, path_prefix, scale=False):
    accs = []
    pres = []
    recs = []
    F1s = []
    for l in range(num_bins):
        cm, acc, pre, rec, F1 = cf.confusion_matrix(actual[l], predicted[l])
        accs.append(acc)
        pres.append(pre)
        recs.append(rec)
        F1s.append(F1)
           
    # plot accuracy
    plt.figure()
    plt.xlabel(variable)
    plt.ylabel('Accuracy')
    if scale:
        plt.plot([b for b in bins], accs, 'ko', rasterized=True)
    else:
        plt.plot([str(round(b, 2)) for b in bins], accs, 'ko', rasterized=True)
    filename = '{}_accuracy.pdf'.format(path_prefix, variable)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    
    #plot precision
    plt.figure()
    plt.xlabel(variable)
    plt.ylabel('Precision')
    if scale:
        plt.plot([b for b in bins], pres, 'ko', rasterized=True)
    else:
        plt.plot([str(round(b, 2)) for b in bins], pres, 'ko', rasterized=True)
    filename = '{}_precision.pdf'.format(path_prefix, variable)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    
    # plot recall
    plt.figure()
    plt.xlabel(variable)
    plt.ylabel('Recall')
    if scale:
        plt.plot([b for b in bins], recs, 'ko', rasterized=True)
    else:
        plt.plot([str(round(b, 2)) for b in bins], recs, 'ko', rasterized=True)
    filename = '{}_recall.pdf'.format(path_prefix, variable)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    
    # plot F1
    plt.figure()
    plt.xlabel(variable)
    plt.ylabel('F1')
    if scale:
        plt.plot([b for b in bins], F1s, 'ko', rasterized=True)
    else:
        plt.plot([str(round(b, 2)) for b in bins], F1s, 'ko', rasterized=True)
    filename = '{}_F1.pdf'.format(path_prefix, variable)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def noise_test(result_foldername):
    if not os.path.exists(result_foldername):
        os.mkdir(result_foldername)
    noise_values = np.array([0,50,100,150,200,300,400,500,750,1000]) 
    altered_noise = noise_values * 10**(-6) + 10**(-9)
    actuals = []
    predicteds = []
    for noise in altered_noise:
        actual, predicted = mf_pipeline(None, result_foldername, noise=noise, mock=True, num_simulations=1000)
        actuals.append(actual)
        predicteds.append(predicted)
    
    plot_completeness('Noise [ppm]', noise_values, actuals, predicteds, 10, "./{}/noise".format(result_foldername), scale=True)
    return actuals, predicteds

def alpha_test(result_foldername):
    if not os.path.exists(result_foldername):
        os.mkdir(result_foldername)
    alpha_values = np.array([.005*i for i in range(1, 11)]) 
    actuals = []
    predicteds = []
    for alpha in alpha_values:
        actual, predicted = mf_pipeline(None, result_foldername, mock=True, num_simulations=1000, threshold=alpha)
        actuals.append(actual)
        predicteds.append(predicted)
    
    plot_completeness('Alpha', alpha_values, actuals, predicteds, 10, "./{}/alpha".format(result_foldername), scale=True)
    return actuals, predicteds
    
def plot_detection(detection, window, flag, pdf, best_template, flare_template, flare_time_window, times, lcur, flat_lcur, correlations, centroid, flat_centroid):
    fig, ax = plt.subplots(5, sharex=True, figsize=(6, 10))
    time_window = times[detection:detection+window]
    
    # plot light curve window
    detection_sample = lcur[detection:detection+window]
    ax[0].set_title('{} Detection at t = {}'.format(flag, round(times[detection],2)))
    ax[0].plot(time_window, detection_sample, 'ko', rasterized=True, markersize=2)
    ax[0].set_ylabel('PDCSAP Flux')
    
    # plot flat light curve window
    flat_detection = mlc.offset_and_normalize(flat_lcur[detection:detection+window])
    ax[1].plot(time_window, flat_detection, 'ko', rasterized=True, markersize=2)
    ax[1].plot(time_window, best_template, 'bo', rasterized=True, markersize=1)
    if flare_template is not None:
        ax[1].plot(flare_time_window, flare_template, 'go', rasterized=True, markersize=1)
        ax[1].legend(["Detection", "Gaussian", "Flare"])
    else:
        ax[1].legend(["Detection", "Gaussian"])
    ax[1].set_ylabel('Relative Flux')
    
    # plot correlation window
    correlation_window = correlations[detection-window//2:detection+window-window//2]
    ax[2].plot(times[detection:detection+window], correlation_window, 'ko', rasterized=True, markersize=2)
    ax[2].set_ylabel('Correlation')
    
    # plot centroid
    centroid_window = centroid[detection:detection+window]
    ax[3].plot(time_window, centroid_window, 'ko', rasterized=True, markersize=2)
    ax[3].set_ylabel('Centroid')
    
    # plot flat centroid
    flat_centroid_window = flat_centroid[detection:detection+window]
    ax[4].plot(time_window, flat_centroid_window, 'ko', rasterized=True, markersize=2)
    ax[4].set_ylabel('Rlative Centroid')
    
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()



foldername = input("Input name of new results folder: ")
os.mkdir(foldername)

templates = []
widths = [5*j + 15 for j in range(10)] + [10*j + 70 for j in range(4)]
for width in widths: 
    templates.append(mlc.offset_and_normalize(mlc.generate_template(1, width)))
    
for sector in globals().get("SECTORS", range(1, 15)):
    flags = ['SL', 'edge', 'transit', 'flare', 'highSL', 'ambiguousSL', 'ambiguousFlare', 'centroid']
    results = {flag: set() for flag in flags}
    if not os.path.exists("{}/Sector{}".format(foldername, sector)):
        os.mkdir("{}/Sector{}".format(foldername, sector))
        os.mkdir("{}/Sector{}/SL_Files".format(foldername, sector))
        for flag in flags:
            os.mkdir("{}/Sector{}/{}".format(foldername, sector, flag))
            
    download_file = "tesscurl_sector_{}_lc.sh".format(sector)
    with open(download_file) as curl_file:
        for command in curl_file:
            split_command = command.split(" ")
            if len(split_command) > 1:
                os.system(command)
                lcur_file = split_command[5]
                tickID = lcur_file.split("-")[2]
                
                try:
                    with fits.open(lcur_file, mode="readonly", memmap=False) as hdulist:
                        tess_bjds = hdulist[1].data['TIME']
                        #sap_fluxes = hdulist[1].data['SAP_FLUX']
                        pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
                        centroid = hdulist[1].data['MOM_CENTR1']
                except:
                    print('Could not open file: {}'.format(lcur_file))
                    continue
                
                # Create a light curve object with the file data and run pipeline on it
                lcur_object = Light_Curve(tess_bjds, pdcsap_fluxes, centroid, templates)
                lcur_object.run_pipeline()
                
                if lcur_object.result:
                    results[flag].add(tickID)
                    folder = '{}/Sector{}/{}'.format(foldername, sector, lcur_object.flag)
                    pdf = pdfs.PdfPages('{}/lcur_{}.pdf'.format(folder, tickID))
                         
                    # create plot for full light curve data
                    fig, ax = plt.subplots(5, sharex=True, figsize=(6, 10))
                    #fig.suptitle(filename)
                    
                    # plot light curve
                    ax[0].set_title("Light Curve {}".format(tickID))
                    ax[0].plot(lcur_object.valid_times, lcur_object.valid_fluxes, 'ko', rasterized=True, markersize=1)
                    ax[0].set_ylabel('PDCSAP Flux')
                    
                    # plot flat light curve
                    ax[1].plot(lcur_object.valid_times, lcur_object.flat_lcur, 'ko', rasterized=True, markersize=1)
                    ax[1].set_ylabel('Relative Flux')
                    
                    # plot correlation
                    ax[2].plot(lcur_object.valid_times[:len(lcur_object.correlations)], lcur_object.correlations, 'ko', rasterized=True, markersize=1)
                    ax[2].plot([lcur_object.valid_times[0], lcur_object.valid_times[len(lcur_object.correlations)-1]], [Light_Curve.alpha, Light_Curve.alpha], '--', color='orange', rasterized=True)
                    if lcur_object.flag == 'highSL':
                        ax[2].plot([lcur_object.valid_times[0], lcur_object.valid_times[len(lcur_object.correlations)-1]], [Light_Curve.beta, Light_Curve.beta], 'b--', rasterized=True)
                    ax[2].set_ylabel('Correlation')
                      
                    # plot centroid
                    ax[3].plot(lcur_object.valid_times, lcur_object.valid_centroid, 'ko', rasterized=True, markersize=1)
                    ax[3].set_ylabel('Centroid')
                    
                    # plot flat centroid
                    ax[4].plot(lcur_object.valid_times, lcur_object.flat_centroid, 'ko', rasterized=True, markersize=1)
                    ax[4].set_ylabel('Relative Centroid')
                    ax[4].set_xlabel('Time [days]')
                    
                    plt.tight_layout()
                    pdf.savefig(fig)
                    plt.close()
                    
                    # zoomed in plot on location of each real positive detection
                    for detection in lcur_object.real_detections:
                        plot_detection(detection[0], lcur_object.window, 'SL', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each ambiguous positive detection
                    for detection in lcur_object.ambiguous_real:
                        plot_detection(detection[0], lcur_object.window, 'Ambiguous SL', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each ambiguous flare detection
                    for detection in lcur_object.ambiguous_flare:
                        plot_detection(detection[0], lcur_object.window, 'Ambiguous Flare', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                        
                    # zoomed in plot on location of each flare detection
                    for detection in lcur_object.flare_detections:
                        plot_detection(detection[0], lcur_object.window, 'Flare', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                    # zoomed in plot on location of each centroid detection
                    for detection in lcur_object.centroid_detections:
                        plot_detection(detection[0], lcur_object.window, 'Centroid', pdf, detection[1], 
                                       detection[2], detection[3], lcur_object.valid_times, lcur_object.valid_fluxes, 
                                       lcur_object.flat_lcur, lcur_object.correlations, lcur_object.valid_centroid, 
                                       lcur_object.flat_centroid)
                    
                    pdf.close()
                    
                if lcur_object.flag in ["SL", "highSL"]:
                    os.system("mv {} {}/Sector{}/SL_Files".format(lcur_file, foldername, sector))
                else:
                    os.system("rm -f {}".format(lcur_file))
                    
        # make a pie chart of the distribution of light curves in each bin
        pie_slices = [len(results[flag]) for flag in flags]
        plt.figure()
        plt.title('Distribution of Positive Detections')
        plt.pie(pie_slices, labels=flags)
        plt.savefig("{}/Sector{}/distribution.pdf".format(foldername, sector))
        plt.close()
                