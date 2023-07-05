"""This code is modified from part of the PyCBC-Tutorials and
 generates an illustration of matched filtering method, using event GW150914 as an example."""
# In[0]:

import os
import matplotlib.pyplot as plt
from pycbc.catalog import Merger
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.waveform import get_td_waveform
from pycbc.filter import resample_to_delta_t, highpass, matched_filter, sigma
from pycbc.waveform.utils import taper_timeseries
from bisect import bisect_right


# In[1]:

merger = Merger("GW150914")
strain = merger.strain('H1')
# downsample the data to 2048Hz, since high_f content is not important here
strain = resample_to_delta_t(highpass(strain, 15.0), 1.0/2048)
# the spike at the boundaries: the filter wraps the data to make it cyclic
conditioned = strain.crop(2, 2)  # Remove 2s of data from both sides

# Use the pycbc.psd.welch method to estimate the psd of this time segment
psd = conditioned.psd(4)
psd = interpolate(psd, conditioned.delta_f)
psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),
                                  low_frequency_cutoff=15)

# Here we assume this is a non-spinning equal mass binary, as the high signal-to-noise ratio(SNR)
# of GW150914, this option doesn't have significant loss of measured SNR
hp, _ = get_td_waveform(approximant="SEOBNRv4_opt", mass1=36, mass2=36,
                        delta_t=conditioned.delta_t, f_lower=20)
hp = taper_timeseries(hp, tapermethod='start')  # make the waveform gradually decreace at start
hp.resize(len(conditioned))  # cutoff in the end
# shift the data so that the merger is approximately at the first bin
template = hp.cyclic_time_shift(hp.start_time)
# time stamps are *not* in general affected

# plt.figure(figsize=(6, 4))
# plt.plot(template)
# plt.tight_layout()

# plt.figure(figsize=(6, 4))
# plt.plot(template.sample_times,template)
# plt.tight_layout()


# In[2]:

snr = matched_filter(template, conditioned, psd=psd, low_frequency_cutoff=20)
# remove 4s at both sides for the PSD filtering
# remove additional 4s at the beginning to account for the template length
# (generous for this short template)
# longer signal like BNS would require much more padding at the beginning
snr = snr.crop(4 + 4, 4)
peak = snr.abs_arg_max()
snrp = snr[peak]
time = snr.sample_times[peak]
print(f'We found a signal at {time:.2f}s with SNR {abs(snrp):.2f}')

plt.figure(figsize=(9, 3))
plt.plot(snr.sample_times, abs(snr), lw=1)
plt.xlim(int(snr.start_time), int(snr.end_time))
plt.ylim(0, 20)
plt.xlabel('Time (s)')
plt.ylabel('Signal-to-Noise Ratio')
plt.tight_layout()
# plt.savefig('MatchedFiltering.pdf')


# In[3]:

dt = time - conditioned.start_time
# scale the template so that it would have SNR 1 in this data
aligned = template/sigma(template, psd=psd, low_frequency_cutoff=20.0)
# Scale the template amplitude and phase to the peak value
aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
# Shift the template to the peak time
aligned = aligned.cyclic_time_shift(dt)
aligned.start_time = conditioned.start_time

# whiten and bandpass
white_template = (aligned.to_frequencyseries() / psd**0.5).to_timeseries()
white_template = white_template.highpass_fir(30, 512).lowpass_fir(300, 512)
white_data = (conditioned.to_frequencyseries() / psd**0.5).to_timeseries()
white_data = white_data.highpass_fir(30, 512).lowpass_fir(300, 512)

fig_dir = 'fig_matched_filtering'
os.makedirs(fig_dir, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex='col', figsize=(12, 5))
ax1.set_ylim(0, 20)
ax1.set_ylabel('Signal-to-Noise Ratio')
ax2.plot(white_data.sample_times, white_data, color='steelblue', lw=1, label="Data")
ax2.set_xlim(merger.time-0.8, merger.time+0.8)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Whitened Strain')

for i in range(321):
    shift = -0.8+0.005*i
    white_template.start_time = white_data.start_time + shift
    idx = bisect_right(snr.sample_times, time + shift)

    snr_p = ax1.plot(snr.sample_times[:idx], abs(snr[:idx]), color='darkblue', lw=1)
    white_p = ax2.plot(white_template.sample_times, white_template, color='darkorange', lw=1, label="Template")
    ax2.legend(loc='upper right')
    # plt.tight_layout()
    plt.savefig(f'{fig_dir}/AlignedWaveform_{i:0>3d}.png')
    snr_p[0].remove()
    white_p[0].remove()


# In[3*]:

white_template.start_time = white_data.start_time
plt.figure(figsize=(9, 3))
plt.plot(white_data.sample_times, white_data, label="Data", lw=1)
plt.plot(white_template.sample_times, white_template, label="Template", lw=1)
plt.xlim(merger.time-0.8, merger.time+0.8)
plt.xlabel('Time (s)')
plt.ylabel('Whitened Strain')
plt.legend(loc='upper right')
plt.tight_layout()
# plt.savefig('AlignedWaveform.pdf')
