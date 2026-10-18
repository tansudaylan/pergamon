import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal 
import rahvar

def P_rng(low=1, high=27):
    logP = np.random.uniform(np.log(low),np.log(high))
    P = np.e**logP
    return P

def mbh_rng():
    low = 5
    mbh = np.random.normal(7.8, 0.6)
    if mbh < low:
        mbh = low
    return mbh

def i_rng():
    cosi  = np.random.random()
    i = math.acos(cosi) * np.random.choice([1,-1])
    return i

def get_amplitudes(P, i, M_BH, M_S=1, R_S=1, rho_S=1.41):
    G = 0.000295913010 # au * 1/solar mass * (au/day)^2
    s_ev = .0189 * math.sin(i)**2 * P**(-2) * (1/rho_S) * (1/(1 + M_S/M_BH))                # unitless 
    s_beam = .0028 * math.sin(i) * P**(-1/3) * (M_BH + M_S)**(-2/3) * M_BH                  # unitless
    a = (G * P**2 * (M_BH + M_S) / (4 * np.pi**2))**(1/3)  # au
    s_sl = 7.15*10**(-5) * R_S**(-2) * P**(2/3) * M_BH * (M_BH + M_S)**(1/3)                # unitless
    tau_sl = .075 * math.pi/4 * P**(1/3) * (M_BH + M_S)**(-1/3) * R_S                       # days
    return s_ev, s_beam, s_sl, tau_sl, a

def generate_light_curve(P, i, M_BH, M_S=1, R_S=1, rho_S=1.41, std=.00006):
    '''
    Takes in various input parameters and generates a 1000 dimensional light curve
    based on EV, SL and beam signals.

    P: Orbital Period (days)
    M_BH: Black Hole Mass (solar masses)
    M_S: Stellar Mass (solar masses)
    R_S: Stellar Radius (solar radii)
    rho_S: Stellar Density (g/cm^3)
    i: orbital inclination (angle)

    returns: size 1000 array containing the light curve
    '''
    
    s_ev, s_beam, s_sl, tau_sl, a = get_amplitudes(P, i, M_BH, M_S=1, R_S=1, rho_S=1.41)

    # signal will represent 27 days
    num_bins = 19440    # 2 minute sampling period
    num_days = 27
    bin_size = num_days/num_bins
    signal = np.zeros(num_bins)
    EV = np.zeros(num_bins)
    Beam = np.zeros(num_bins)
    SL = np.zeros(num_bins)
    t = 0
    t0 = 0
    t0s = []
    while t0 < num_days:
        t0s.append(t0)
        t0 += P
        
    for j in range(num_bins):
        ev =  -s_ev * math.cos(4*math.pi*t/P)
        beam = s_beam * math.sin(2*math.pi*t/P)
        closest_t0 = 0
        difference = t 
        for t0 in t0s:
            new_difference = abs(t-t0)
            if new_difference <= difference:
                closest_t0 = t0
                difference = new_difference
            else:
                break
        sl = rahvar.rahvar_main(t*24, closest_t0*24, i, M_BH, M_S, a)
        #sl = s_sl if t%P <= tau_sl else 0
        signal[j] += ev + beam + sl + np.random.normal(0, std) 
        EV[j] = ev 
        Beam[j] = beam 
        SL[j] = sl  
        t = t + bin_size
    
    return signal, EV+1, Beam+1, SL

def generate_flat_signal(std):
    num_bins = 19440
    return np.array([np.random.normal(0, std) for _ in range(num_bins)])

def plot_lc(lc, P, M_BH, i, filename=None, num_days=27, num_bins=19440, EV=None, Beam=None, SL=None):
    bin_size = num_days/num_bins
    plt.figure()
    plt.xlabel('Time [days]')
    plt.ylabel('Relative Flux')
    title = 'P = ' + str(round(P,2)) + ' days, ' r'$M_{BH} = $' + str(round(M_BH,2)) + r' $ M_{\odot}, cosi = $' + str(round(math.cos(i),2))
    plt.title(title)
    lc_plot = plt.plot([i*bin_size for i in range(num_bins)], lc, 'kor', label='Signal', rasterized=True)
    handles = lc_plot
    if EV is not None:
        EV_plot = plt.plot([i*bin_size for i in range(num_bins)], EV, 'b--', label='EV', rasterized=True)
        handles += EV_plot
    if Beam is not None:
        Beam_plot = plt.plot([i*bin_size for i in range(num_bins)], Beam, 'g--', label='Beam', rasterized=True)
        handles += Beam_plot
    if SL is not None:
        SL_plot = plt.plot([i*bin_size for i in range(num_bins)], SL, 'r--', label='SL', rasterized=True)
        handles += SL_plot

    plt.legend(handles=handles, loc="upper right")
    
    if filename is not None:
        plt.tight_layout()
        plt.savefig(filename)
        
    plt.close()
    
def plot_sl(lc, P, M_BH, i, filename=None, num_days=27, num_bins=19440, EV=None, Beam=None, SL=None):
    bin_size = num_days/num_bins
    plt.figure()
    plt.xlabel('Time [days]')
    plt.ylabel('Relative Flux')
    #title = 'P = ' + str(round(P,2)) + ' days, ' r'$M_{BH} = $' + str(round(M_BH,2)) + r' $ M_{\odot}, cosi = $' + str(round(math.cos(i),2))
    #plt.title(title)
    #lc_plot = plt.plot([i*bin_size for i in range(num_bins)], lc, 'k', label='Signal')
    handles = []
# =============================================================================
#     if EV is not None:
#         EV_plot = plt.plot([i*bin_size for i in range(num_bins)], EV, 'b--', label='EV')
#         handles += EV_plot
#     if Beam is not None:
#         Beam_plot = plt.plot([i*bin_size for i in range(num_bins)], Beam, 'g--', label='Beam')
#         handles += Beam_plot
# =============================================================================
    if SL is not None:
        SL_plot = plt.plot([i*bin_size for i in range(num_bins)], SL, 'r--', label='SL')
        handles += SL_plot

    plt.legend(handles=handles, loc="upper right")
    
    if filename is not None:
        plt.tight_layout()
        plt.savefig(filename)
        
    plt.close()
    
def supersample(signal):
    '''
    Supersample the signal by averaging out each data point with surrounding data points. 
    Uses a 15 sample window to eliminate noise.
    
    signal: np array, signal to be supersampled
    
    returns: np array, supersampled signal 
    '''
    
    # supersample with a 30 min cadence 
    num_bins = len(signal)
    total = 0
    for k in range(15):
        total += signal[k]
        
    for l in range(num_bins-8):
        if  l > 6:
            signal[l] = total/15
            total += signal[l+8] - signal[l-7]
            
    return signal


def upsample(signal, up):
    '''
    Takes in a signal and upsamples it by an upsampling factor up
    
    signal: numpy array or list, signal to be upsampled
    up: int (>1), upsampling factor
    
    returns: np array, upsampled signal
    '''
    
    if up < 1 or type(up) != int:
        raise TypeError("upsampling factor not an integer greater than 1")
        
    new_signal = []
        
    for i in range(len(signal)-1):
        new_signal.append(signal[i])
        dy = (signal[i+1] - signal[i])/up
        for j in range(1, up):
            new_signal.append(signal[i] + j*dy)
    
    last_dy = (signal[-1] - signal[-2])/up
    new_signal.append(signal[-1])
    
    for j in range(1, up):
        new_signal.append(signal[-1] + j*last_dy)
    
    return np.array(new_signal)


def downsample(signal, down):
    '''
    Takes in a signal and downsamples it by a downsampling factor down
    
    signal: numpy array or list, signal to be upsampled
    down: int (>1), downsampling factor
    
    returns: np array, downsampled signal
    '''
    
    if down<1 or type(down) != int:
        raise TypeError("upsampling factor not an integer greater than 1")
    
    new_signal = []
    
    for i in range(len(signal)):
        if i%down == 0:
            new_signal.append(signal[i])
            
    return new_signal

def resample(signal, desired_length):
    '''
    Takes in a signal and resamples it to the desired length by upsampling and downsampling
    
    signal: numpy array or list, signal to be resampled
    desired_length: int, length of the resampled signal
    
    returns: np array, resampled signal
    '''
    
    signal_length = len(signal)
    LCM = signal_length * desired_length/ math.gcd(signal_length, desired_length)
    up_factor = int(LCM/signal_length)
    down_factor = int(LCM/desired_length)
    up_signal = upsample(signal, up_factor)
    result = downsample(up_signal, down_factor)
    return result


def offset_and_normalize(signal):
    '''
    Takes a signal and normalizes it by subtracting the mean and dividing by the magnitude
    
    signal: numpy array or list, signal to be resampled
    
    returns: np array, normalized signal
    '''
    
    avg = np.mean(signal)
    total = 0
    
    for num in signal:
        total += (num - avg)**2
        
    mag = total**(.5)
    
    return (signal - avg) / mag


def correlation(signal_1, signal_2):
    '''
    Takes two signals and finds the correlation between them 
    
    signal_1, signal_2: np arrays, signals to be correlated
    
    returns: int (0 to 1), correlation between the two signals
    '''
    
# =============================================================================
#     norm_signal_1 = offset_and_normalize(signal_1)
#     norm_signal_2 = offset_and_normalize(signal_2)
#     
#     return np.sum(norm_signal_1 * norm_signal_2)
# =============================================================================
    return np.sum(signal_1 * signal_2)
    

def flatten(lc, P, i, M_BH, M_S=1, R_S=1, rho_S=1.41, num_days=27):
    '''
    Takes a light curve and flattens it out by subtracting out the ellipsoidal variation and
    beaming signals. Flattened curve contains self lensing aplicfication only
    
    lc: np array, light curve to be flattened
    P: int or float, Orbital period (days)
    s_ev: float, ellipsoidal variation amplitude
    s_beam: float, beaming amplitude
    num_days: int or float, length of signal in days
    
    returns: np array, flattened light curve
    '''
    
    s_ev, s_beam, s_sl, tau_sl = get_amplitudes(P, i, M_BH, M_S, R_S, rho_S)
    flattened_lc = lc.copy()
    num_bins = len(lc)    # 2 minute sampling period
    bin_size = num_days/num_bins
    t = 0
    
    for j in range(num_bins):
        ev =  -s_ev * math.cos(4*math.pi*t/P)
        beam = s_beam * math.sin(2*math.pi*t/P)
        flattened_lc[j] -= ev + beam
        t += bin_size
        
    return flattened_lc
    
def match_filter(lc, template, P=None, i=None, M_BH=None, M_S=1, R_S=1, rho_S=1.41, num_days=27, threshold=.5, mock=False, alpha=.5):
    '''
    Runs a template through a signal and correlates the two to find if the two signals have a correlation 
    that exceeds a given threshold.
    
    lc: np array, light curve
    template: np array, template to be correlated (None if mock)
    P: Orbital Period (days)
    i: Orbital Inclination (radians)
    M_BH: Black Hole Mass (solar masses)
    threshold: float, value which correlation much exceed for a positive classification
    mock: bool, whether or not the signal is mock data
    
    returns: bool and float, whether or not the signal matches the template and correlation
    between the signal and the template in the given window
    '''

    num_bins = len(lc)
    
    if mock:
        s_ev, s_beam, s_sl, tau_sl, a = get_amplitudes(P, i, M_BH, M_S, R_S, rho_S)
        bin_size = num_days/num_bins
        sl_bins = int(tau_sl // bin_size)
        template = generate_template(s_sl, sl_bins)
        threshold = s_sl**2#*sl_bins*alpha
                
    window = len(template)
    highest_corr = 0
    correlations = np.zeros(num_bins - window)
    num_sl_events = 0
    last_corr = False
    for j in range(num_bins-window):
        corr = correlation(lc[j:j+window], template)
        correlations[j] = corr
        
        if corr > threshold:
            if not last_corr:
                num_sl_events += 1
            
        last_corr = corr > threshold
        
        if corr > highest_corr:
            highest_corr = corr
        
    return {"result": highest_corr > threshold, "highest_corr": highest_corr, "correlations": correlations, 
            "template": template, "threshold": threshold, "sl_events": num_sl_events}

def plot_corr(correlations, P, M_BH, i, alpha, window, threshold, filename, num_days=27, num_bins=19440):
    plt.figure()
    plt.xlabel('Time [days]')
    plt.ylabel('Correlation')
    title = 'P = ' + str(round(P,2)) + ' [days], ' r'$M_{BH} = $' + str(round(M_BH,2)) + r' $ [M_{\odot}], cosi = $' + str(round(math.cos(i),2)) + r', $ \alpha = $' + str(round(alpha, 2))
    plt.title(title)
    plt.plot([x*num_days/num_bins for x in range(num_bins-window)], correlations, 'k', rasterized=True)
    plt.plot([x*num_days/num_bins for x in range(num_bins-window)], [threshold for _ in range(num_bins-window)], 'b--', rasterized=True)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def generate_template(s_sl, sl_bins):
    '''
    Creates a square pulse with the correct length to match the predicted sl period 
    of a light curve.  Square pulse is meant to match sl amplitude.
    
    s_sl: float, s_sl amplitude
    sl_bins: int, length of 
    '''
    def gaussian(x, b):
        w = sl_bins
        a = s_sl
        return a * math.e**(-(x-b)**2 / (w**2))
        
    return np.array([gaussian(3*x, 4.5*sl_bins) for x in range(3*sl_bins)])
    