// basic_fir_design.c : design basic FIR filter

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

int main() {
    // parameters and simulation options
    unsigned int h_len =   57;  // filter length (samples)
    float        fc    = 0.20;  // normalized cutoff frequency
    unsigned int nfft  =  800;  // 'FFT' size (actually DFT)

    unsigned int i;
    unsigned int k;

    // design filter
    float h[h_len];
    for (i=0; i < h_len; i++) {
        // generate time vector, centered at zero
        float t = (float)i + 0.5f - 0.5f*(float)h_len;

        // generate sinc function (time offset in 't' prevents divide by zero)
        float s = sinf(2*M_PI*fc*t + 1e-6f) / (2*M_PI*fc*t + 1e-6f);

        // generate Hamming window
        float w = 0.53836 - 0.46164*cosf((2*M_PI*(float)i)/((float)(h_len-1)));

        // generate composite filter coefficient
        h[i] = s * w;
    }

    // print line legend to standard output
    printf("# %12s %12s\n", "frequency", "PSD [dB]");

    // run filter analysis with discrete Fourier transform
    for (i=0; i < nfft; i++) {
        // accumulate energy in frequency (also apply frequency shift)
        float complex H = 0.0f;
        float frequency = (float)i/(float)nfft-0.5f;
        for (k=0; k < h_len; k++)
            H += h[k] * cexpf(_Complex_I*2*M_PI*k*frequency);

        // print resulting power spectral density to standard output
        printf("  %12.8f %12.6f\n", frequency, 20*log10(cabsf(H*2*fc)));
    }

    return 0;
}