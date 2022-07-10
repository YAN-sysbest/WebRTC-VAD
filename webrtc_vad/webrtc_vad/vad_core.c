/*
 *  Copyright (c) 2011 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/*
 * This file includes the implementation of the core functionality in VAD.
 * For function description, see vad_core.h.
 */

#include <assert.h>
#include "vad_core.h"

#include "signal_processing_library.h"//差三个inline函数
#include "typedefs.h"
#include "vad_defines.h"
//#include "vad_filterbank.h"
//#include "vad_gmm.h"
//#include "vad_sp.h"

// Spectrum Weighting
static const WebRtc_Word16 kSpectrumWeight[6] = { 6, 8, 10, 12, 14, 16 };
static const WebRtc_Word16 kNoiseUpdateConst = 655; // Q15
static const WebRtc_Word16 kSpeechUpdateConst = 6554; // Q15
static const WebRtc_Word16 kBackEta = 154; // Q8
// Minimum difference between the two models, Q5
static const WebRtc_Word16 kMinimumDifference[6] = {
    544, 544, 576, 576, 576, 576 };
// Upper limit of mean value for speech model, Q7
static const WebRtc_Word16 kMaximumSpeech[6] = {
    11392, 11392, 11520, 11520, 11520, 11520 };
// Minimum value for mean value
static const WebRtc_Word16 kMinimumMean[2] = { 640, 768 };
// Upper limit of mean value for noise model, Q7
static const WebRtc_Word16 kMaximumNoise[6] = {
    9216, 9088, 8960, 8832, 8704, 8576 };
// Start values for the Gaussian models, Q7
// Weights for the two Gaussians for the six channels (noise)
static const WebRtc_Word16 kNoiseDataWeights[12] = {
    34, 62, 72, 66, 53, 25, 94, 66, 56, 62, 75, 103 };
// Weights for the two Gaussians for the six channels (speech)
static const WebRtc_Word16 kSpeechDataWeights[12] = {
    48, 82, 45, 87, 50, 47, 80, 46, 83, 41, 78, 81 };
// Means for the two Gaussians for the six channels (noise)
static const WebRtc_Word16 kNoiseDataMeans[12] = {
    6738, 4892, 7065, 6715, 6771, 3369, 7646, 3863, 7820, 7266, 5020, 4362 };
// Means for the two Gaussians for the six channels (speech)
static const WebRtc_Word16 kSpeechDataMeans[12] = {
    8306, 10085, 10078, 11823, 11843, 6309, 9473, 9571, 10879, 7581, 8180, 7483
};
// Stds for the two Gaussians for the six channels (noise)
static const WebRtc_Word16 kNoiseDataStds[12] = {
    378, 1064, 493, 582, 688, 593, 474, 697, 475, 688, 421, 455 };
// Stds for the two Gaussians for the six channels (speech)
static const WebRtc_Word16 kSpeechDataStds[12] = {
    555, 505, 567, 524, 585, 1231, 509, 828, 492, 1540, 1079, 850 };

static const int kInitCheck = 42;

//vad_gmm

static const int32_t kCompVar = 22005;
static const int16_t kLog2Exp = 5909;  // log2(exp(1)) in Q12.

int32_t WebRtcSpl_DivW32W16(int32_t num, int16_t den)
{
    if (den != 0)
    {
        return (int32_t)(num / den);
    }
    else
    {
        return (int32_t)0x7FFFFFFF;
    }
}

int32_t WebRtcVad_GaussianProbability(int16_t input,
    int16_t mean,
    int16_t std,
    int16_t* delta) {
    int16_t tmp16, inv_std, inv_std2, exp_value = 0;
    int32_t tmp32;

    tmp32 = (int32_t)131072 + (int32_t)(std >> 1);
    inv_std = (int16_t)WebRtcSpl_DivW32W16(tmp32, std);

    tmp16 = (inv_std >> 2);  // Q10 -> Q8.
    inv_std2 = (int16_t)((int32_t)(tmp16 * tmp16) >> 2);

    tmp16 = (input << 3);  // Q4 -> Q7
    tmp16 = tmp16 - mean;  // Q7 - Q7 = Q7

    *delta = (int16_t)((int32_t)(inv_std2 * tmp16) >> 10);

    tmp32 = ((int32_t)((*delta) * tmp16) >> 9);
    if (tmp32 < 22005) {
        tmp16 = (int16_t)((int32_t)(5909 * (int16_t)tmp32) >> 12);

        tmp16 = -tmp16;
        exp_value = (0x0400 | (tmp16 & 0x03FF));
        tmp16 ^= 0xFFFF;
        tmp16 >>= 10;
        tmp16 += 1;
        exp_value >>= tmp16;
    }
    return ((int32_t)(((int16_t)(inv_std)) * ((int16_t)(exp_value))));
}


//vad_sp

static const int16_t kAllPassCoefsQ13[2] = { 5243, 1392 };  // Q13
#define WEBRTC_SPL_WORD16_MAX       32767
#define WEBRTC_SPL_WORD16_MIN       -32768
 
void WebRtcVad_Downsampling(int16_t* signal_in,
    int16_t* signal_out,
    int32_t* filter_state,
    int in_length) {
    int16_t tmp16_1 = 0, tmp16_2 = 0;
    int32_t tmp32_1 = filter_state[0];
    int32_t tmp32_2 = filter_state[1];
    int n = 0;
    int half_length = (in_length >> 1);  // Downsampling by 2 gives half length.

    // Filter coefficients in Q13, filter state in Q0.
    for (n = 0; n < half_length; n++) {
        // All-pass filtering upper branch.
        //tmp16_1 = (int16_t) ((tmp32_1 >> 1) +  WEBRTC_SPL_MUL_16_16_RSFT(kAllPassCoefsQ13[0], *signal_in, 14));
        tmp16_1 = (int16_t)((tmp32_1 >> 1) + ((int32_t)kAllPassCoefsQ13[0] * (*signal_in) >> 14));
        *signal_out = tmp16_1;
        //tmp32_1 = (int32_t) (*signal_in++) -        WEBRTC_SPL_MUL_16_16_RSFT(kAllPassCoefsQ13[0], tmp16_1, 12);
        tmp32_1 = (int32_t)(*signal_in++) - ((int32_t)kAllPassCoefsQ13[0] * tmp16_1 >> 12);
        // All-pass filtering lower branch.
        //tmp16_2 = (int16_t) ((tmp32_2 >> 1) +        WEBRTC_SPL_MUL_16_16_RSFT(kAllPassCoefsQ13[1], *signal_in, 14));
        tmp16_2 = (int16_t)((tmp32_2 >> 1) + ((int32_t)kAllPassCoefsQ13[1] * (*signal_in) >> 14));
        *signal_out++ += tmp16_2;
        tmp32_2 = (int32_t)(*signal_in++) - ((int32_t)kAllPassCoefsQ13[1] * tmp16_2 >> 12);
    }
    // Store the filter states.
    filter_state[0] = tmp32_1;
    filter_state[1] = tmp32_2;
}
 
int16_t WebRtcVad_FindMinimum(VadInstT* self,
    int16_t feature_value,
    int channel) {
    int i = 0, j = 0;
    int position = -1;
    // Offset to beginning of the 16 minimum values in memory.
    int offset = (channel << 4);
    int16_t current_median = 1600;
    int16_t alpha = 0;
    int32_t tmp32 = 0;
    // Pointer to memory for the 16 minimum values and the age of each value of
    // the |channel|.
    int16_t* age_ptr = &self->index_vector[offset];
    int16_t* value_ptr = &self->low_value_vector[offset];
    int16_t* p1, * p2, * p3;

    assert(channel < NUM_CHANNELS);

    // Each value in |low_value_vector| is getting 1 loop older.
    // Update age of each value in |age_ptr|, and remove old values.
    for (i = 0; i < 16; i++) {
        p3 = age_ptr + i;
        if (*p3 != 100) {
            *p3 += 1;
        }
        else {
            p1 = value_ptr + i + 1;
            p2 = p3 + 1;
            for (j = i; j < 16; j++) {
                *(value_ptr + j) = *p1++;
                *(age_ptr + j) = *p2++;
            }
            *(age_ptr + 15) = 101;
            *(value_ptr + 15) = 10000;
        }
    }

    // Check if |feature_value| is smaller than any of the values in
    // |low_value_vector|. If so, find the |position| where to insert the new
    // value.
    if (feature_value < *(value_ptr + 7)) {
        if (feature_value < *(value_ptr + 3)) {
            if (feature_value < *(value_ptr + 1)) {
                if (feature_value < *value_ptr) {
                    position = 0;
                }
                else {
                    position = 1;
                }
            }
            else if (feature_value < *(value_ptr + 2)) {
                position = 2;
            }
            else {
                position = 3;
            }
        }
        else if (feature_value < *(value_ptr + 5)) {
            if (feature_value < *(value_ptr + 4)) {
                position = 4;
            }
            else {
                position = 5;
            }
        }
        else if (feature_value < *(value_ptr + 6)) {
            position = 6;
        }
        else {
            position = 7;
        }
    }
    else if (feature_value < *(value_ptr + 15)) {
        if (feature_value < *(value_ptr + 11)) {
            if (feature_value < *(value_ptr + 9)) {
                if (feature_value < *(value_ptr + 8)) {
                    position = 8;
                }
                else {
                    position = 9;
                }
            }
            else if (feature_value < *(value_ptr + 10)) {
                position = 10;
            }
            else {
                position = 11;
            }
        }
        else if (feature_value < *(value_ptr + 13)) {
            if (feature_value < *(value_ptr + 12)) {
                position = 12;
            }
            else {
                position = 13;
            }
        }
        else if (feature_value < *(value_ptr + 14)) {
            position = 14;
        }
        else {
            position = 15;
        }
    }

    // If we have a new small value, put it in the correct position and shift
    // larger values up.
    if (position > -1) {
        for (i = 15; i > position; i--) {
            j = i - 1;
            *(value_ptr + i) = *(value_ptr + j);
            *(age_ptr + i) = *(age_ptr + j);
        }
        *(value_ptr + position) = feature_value;
        *(age_ptr + position) = 1;
    }

    // Get |current_median|.
    if (self->frame_counter > 2) {
        current_median = *(value_ptr + 2);
    }
    else if (self->frame_counter > 0) {
        current_median = *value_ptr;
    }

    // Smooth the median value.
    if (self->frame_counter > 0) {
        if (current_median < self->mean_value[channel]) {
            alpha = (int16_t)ALPHA1;  // 0.2 in Q15.
        }
        else {
            alpha = (int16_t)ALPHA2;  // 0.99 in Q15.
        }
    }

    // tmp32 = WEBRTC_SPL_MUL_16_16(alpha + 1, self->mean_value[channel]);
     //tmp32 += WEBRTC_SPL_MUL_16_16(WEBRTC_SPL_WORD16_MAX - alpha, current_median);

    tmp32 = ((int32_t)(((int16_t)(alpha + 1)) * ((int16_t)(self->mean_value[channel]))));
    tmp32 += ((int32_t)(((int16_t)(WEBRTC_SPL_WORD16_MAX - alpha)) * ((int16_t)(current_median))));

    tmp32 += 16384;
    self->mean_value[channel] = (int16_t)(tmp32 >> 15);

    return self->mean_value[channel];
}



//vad_filter

// Constant 160*log10(2) in Q9
static const int16_t kLogConst = 24660;

// Coefficients used by WebRtcVad_HpOutput, Q14
static const int16_t kHpZeroCoefs[3] = { 6631, -13262, 6631 };
static const int16_t kHpPoleCoefs[3] = { 16384, -7756, 5620 };

// Allpass filter coefficients, upper and lower, in Q15
// Upper: 0.64, Lower: 0.17
static const int16_t kAllPassCoefsQ15[2] = { 20972, 5571 };

// Adjustment for division with two in WebRtcVad_SplitFilter
static const int16_t kOffsetVector[6] = { 368, 368, 272, 176, 176, 176 };
/*
void WebRtcVad_HpOutput(int16_t* in_vector,
                        int in_vector_length,
                        int16_t* filter_state,
                        int16_t* out_vector) {
  int i;
  int16_t* in_ptr = in_vector;
  int16_t* out_ptr = out_vector;
  int32_t tmp32 = 0;
  for (i = 0; i < in_vector_length; i++) {
    // all-zero section (filter coefficients in Q14)
    tmp32 = (int32_t) WEBRTC_SPL_MUL_16_16(kHpZeroCoefs[0], (*in_ptr));
    tmp32 += (int32_t) WEBRTC_SPL_MUL_16_16(kHpZeroCoefs[1], filter_state[0]);
    tmp32 += (int32_t) WEBRTC_SPL_MUL_16_16(kHpZeroCoefs[2],
                                            filter_state[1]);  // Q14
    filter_state[1] = filter_state[0];
    filter_state[0] = *in_ptr++;

    // all-pole section
    tmp32 -= (int32_t) WEBRTC_SPL_MUL_16_16(kHpPoleCoefs[1],
                                            filter_state[2]);  // Q14
    tmp32 -= (int32_t) WEBRTC_SPL_MUL_16_16(kHpPoleCoefs[2], filter_state[3]);
    filter_state[3] = filter_state[2];
    filter_state[2] = (int16_t) WEBRTC_SPL_RSHIFT_W32 (tmp32, 14);
    *out_ptr++ = filter_state[2];
  }
}
*/

void WebRtcVad_HpOutput(int16_t* in_vector,
    int in_vector_length,
    int16_t* filter_state,
    int16_t* out_vector) {
    int i;
    int16_t* in_ptr = in_vector;
    int16_t* out_ptr = out_vector;
    int32_t tmp32 = 0;

    for (i = 0; i < in_vector_length; i++) {
        // all-zero section (filter coefficients in Q14)
        tmp32 = (int32_t)(kHpZeroCoefs[0] * (*in_ptr));
        tmp32 += (int32_t)(kHpZeroCoefs[1] * filter_state[0]);
        tmp32 += (int32_t)(kHpZeroCoefs[2] * filter_state[1]);  // Q14
        filter_state[1] = filter_state[0];
        filter_state[0] = *in_ptr++;

        // all-pole section
        tmp32 -= (int32_t)(kHpPoleCoefs[1] * filter_state[2]);  // Q14
        tmp32 -= (int32_t)(kHpPoleCoefs[2] * filter_state[3]);
        filter_state[3] = filter_state[2];
        filter_state[2] = (int16_t)(tmp32 >> 14);
        *out_ptr++ = filter_state[2];
    }
}
/*

void WebRtcVad_Allpass(int16_t* in_vector,
                       int16_t filter_coefficients,
                       int vector_length,
                       int16_t* filter_state,
                       int16_t* out_vector) {

  int i;
  int16_t tmp16 = 0;
  int32_t tmp32 = 0, in32 = 0;
  int32_t state32 = WEBRTC_SPL_LSHIFT_W32((int32_t) (*filter_state), 16); // Q31

  for (i = 0; i < vector_length; i++) {
    tmp32 = state32 + WEBRTC_SPL_MUL_16_16(filter_coefficients, (*in_vector));
    tmp16 = (int16_t) WEBRTC_SPL_RSHIFT_W32(tmp32, 16);
    *out_vector++ = tmp16;
    in32 = WEBRTC_SPL_LSHIFT_W32(((int32_t) (*in_vector)), 14);
    state32 = in32 - WEBRTC_SPL_MUL_16_16(filter_coefficients, tmp16);
    state32 = WEBRTC_SPL_LSHIFT_W32(state32, 1);
    in_vector += 2;
  }

  *filter_state = (int16_t) WEBRTC_SPL_RSHIFT_W32(state32, 16);
}
*/

void WebRtcVad_Allpass(int16_t* in_vector,
    int16_t filter_coefficients,
    int vector_length,
    int16_t* filter_state,
    int16_t* out_vector) {
    int i;
    int32_t tmp32 = 0, in32 = 0;
    int32_t state32 = ((int32_t)((*filter_state) << 16)); // Q31

    for (i = 0; i < vector_length; i++) {
        tmp32 = state32 + ((int32_t)(((int16_t)(filter_coefficients)) * ((int16_t)(*in_vector))));

        *out_vector++ = (int16_t)(tmp32 >> 16);
        state32 = ((int32_t)((*in_vector) << 14)) -
            ((int32_t)(((int16_t)(filter_coefficients)) * ((int16_t)(tmp32 >> 16))));
        state32 = (state32 << 1);
        in_vector += 2;
    }

    *filter_state = (int16_t)(state32 >> 16);
}





void WebRtcVad_SplitFilter(int16_t* in_vector,
    int in_vector_length,
    int16_t* upper_state,
    int16_t* lower_state,
    int16_t* out_vector_hp,
    int16_t* out_vector_lp) {
    int16_t tmp_out;
    int i;
    int half_length = in_vector_length >> 1;

    // All-pass filtering upper branch
    WebRtcVad_Allpass(&in_vector[0], kAllPassCoefsQ15[0], half_length,
        upper_state, out_vector_hp);

    // All-pass filtering lower branch
    WebRtcVad_Allpass(&in_vector[1], kAllPassCoefsQ15[1], half_length,
        lower_state, out_vector_lp);

    // Make LP and HP signals
    for (i = 0; i < half_length; i++) {
        tmp_out = *out_vector_hp;
        *out_vector_hp++ -= *out_vector_lp;
        *out_vector_lp++ += tmp_out;
    }
}
/*

void WebRtcVad_SplitFilter(int16_t* in_vector,
                           int in_vector_length,
                           int16_t* upper_state,
                           int16_t* lower_state,
                           int16_t* out_vector_hp,
                           int16_t* out_vector_lp) {
  int16_t tmp_out;
  int i;
  int half_length = WEBRTC_SPL_RSHIFT_W16(in_vector_length, 1);

  // All-pass filtering upper branch
  WebRtcVad_Allpass(&in_vector[0], kAllPassCoefsQ15[0], half_length,
                    upper_state, out_vector_hp);

  // All-pass filtering lower branch
  WebRtcVad_Allpass(&in_vector[1], kAllPassCoefsQ15[1], half_length,
                    lower_state, out_vector_lp);

  // Make LP and HP signals
  for (i = 0; i < half_length; i++) {
    tmp_out = *out_vector_hp;
    *out_vector_hp++ -= *out_vector_lp;
    *out_vector_lp++ += tmp_out;
  }
}
*/


//get_scaling.c

int WebRtcSpl_GetScalingSquare(int16_t* in_vector, int in_vector_length,
    int times)
{
    int nbits = WebRtcSpl_GetSizeInBits(times);
    int i;
    int16_t smax = -1;
    int16_t sabs;
    int16_t* sptr = in_vector;
    int t;
    int looptimes = in_vector_length;

    for (i = looptimes; i > 0; i--)
    {
        sabs = (*sptr > 0 ? *sptr++ : -*sptr++);
        smax = (sabs > smax ? sabs : smax);
    }
    t = WebRtcSpl_NormW32((int32_t)smax * smax);

    if (smax == 0) return 0;
    else return (t > nbits) ? 0 : nbits - t;
}


//energy.c

int32_t WebRtcSpl_Energy(int16_t* vector, int vector_length, int* scale_factor)
{
    int32_t en = 0;
    int i;
    int scaling = WebRtcSpl_GetScalingSquare(vector, vector_length, vector_length);
    int looptimes = vector_length;
    int16_t* vectorptr = vector;

    for (i = 0; i < looptimes; i++)
    {
        // en += WEBRTC_SPL_MUL_16_16_RSFT(*vectorptr, *vectorptr, scaling);
        en += (int32_t)(*vectorptr) * (*vectorptr) >> scaling;
        vectorptr++;
    }
    *scale_factor = scaling;

    return en;
}


void WebRtcVad_LogOfEnergy(int16_t* vector,
    int vector_length,
    int16_t offset,
    int16_t* power,
    int16_t* log_energy) {
    int shfts = 0, shfts2 = 0;
    int16_t energy_s16 = 0;
    int16_t zeros = 0, frac = 0, log2 = 0;
    int32_t energy = WebRtcSpl_Energy(vector, vector_length, &shfts);

    if (energy > 0) {

        shfts2 = 16 - WebRtcSpl_NormW32(energy);
        shfts += shfts2;
        //energy_s16 = (int16_t) WEBRTC_SPL_SHIFT_W32(energy, -shfts2);
        energy_s16 = (int16_t)(((-shfts2) >= 0) ? ((energy) << (-shfts2)) : ((energy) >> (shfts2)));

        // #define WEBRTC_SPL_SHIFT_W32(x, c) \
         (((c) >= 0) ? ((x) << (c)) : ((x) >> (-(c)))) 


        zeros = WebRtcSpl_NormU32(energy_s16);
        frac = (int16_t)(((uint32_t)((int32_t)(energy_s16) << zeros)
            & 0x7FFFFFFF) >> 21);
        log2 = (int16_t)(((31 - zeros) << 10) + frac);
        //#define WEBRTC_SPL_MUL_16_16_RSFT(a, b, c) \
        (WEBRTC_SPL_MUL_16_16(a, b) >> (c))
      //  *log_energy = (int16_t) WEBRTC_SPL_MUL_16_16_RSFT(kLogConst, log2, 19) + (int16_t) WEBRTC_SPL_MUL_16_16_RSFT(shfts, kLogConst, 9);
        *log_energy = (int16_t)((int32_t)(kLogConst * log2) >> 19) + (int16_t)((int32_t)(shfts * kLogConst) >> 9);


        if (*log_energy < 0) {
            *log_energy = 0;
        }
    }
    else {
        *log_energy = 0;
        shfts = -15;
        energy_s16 = 0;
    }

    *log_energy += offset;

    // Total power in frame,MIN_ENERGY=10
    if (*power <= 10) {
        if (shfts > 0) {
            *power += 10 + 1;
            //#define WEBRTC_SPL_SHIFT_W16(x, c) \
          (((c) >= 0) ? ((x) << (c)) : ((x) >> (-(c))))
        } //else if (WEBRTC_SPL_SHIFT_W16(energy_s16, shfts) > 10) {(energy_s16, shfts) 
        else if ((shfts >= 0) ? ((energy_s16) << shfts) : (energy_s16 >> (-(shfts))) > 10) {
            *power += 10 + 1;
        }
        else {
            //*power += WEBRTC_SPL_SHIFT_W16(energy_s16, shfts);
            *power += (shfts >= 0) ? ((energy_s16) << shfts) : (energy_s16 >> (-(shfts)));
        }
    }
}

//vad_filter主要的函数
int16_t WebRtcVad_get_features(VadInstT* inst,
    int16_t* in_vector,
    int frame_size,
    int16_t* out_vector) {
    int16_t power = 0;
    int16_t hp_120[120], lp_120[120];
    int16_t hp_60[60], lp_60[60];
    // Initialize variables for the first SplitFilter().
    int length = frame_size;
    int frequency_band = 0;
    int16_t* in_ptr = in_vector;
    int16_t* hp_out_ptr = hp_120;
    int16_t* lp_out_ptr = lp_120;

    // Split at 2000 Hz and downsample
    WebRtcVad_SplitFilter(in_ptr, length, &inst->upper_state[frequency_band],
        &inst->lower_state[frequency_band], hp_out_ptr,
        lp_out_ptr);

    // Split at 3000 Hz and downsample
    frequency_band = 1;
    in_ptr = hp_120;
    hp_out_ptr = hp_60;
    lp_out_ptr = lp_60;
    //  ((x) >> (c));length = WEBRTC_SPL_RSHIFT_W16(frame_size, 1);
    length = (frame_size >> 1);
    WebRtcVad_SplitFilter(in_ptr, length, &inst->upper_state[frequency_band],
        &inst->lower_state[frequency_band], hp_out_ptr,
        lp_out_ptr);

    // Energy in 3000 Hz - 4000 Hz
    length = (length >> 1);
    WebRtcVad_LogOfEnergy(hp_60, length, kOffsetVector[5], &power,
        &out_vector[5]);

    // Energy in 2000 Hz - 3000 Hz
    WebRtcVad_LogOfEnergy(lp_60, length, kOffsetVector[4], &power,
        &out_vector[4]);

    // Split at 1000 Hz and downsample
    frequency_band = 2;
    in_ptr = lp_120;
    hp_out_ptr = hp_60;
    lp_out_ptr = lp_60;
    length = (frame_size >> 1);
    WebRtcVad_SplitFilter(in_ptr, length, &inst->upper_state[frequency_band],
        &inst->lower_state[frequency_band], hp_out_ptr,
        lp_out_ptr);

    // Energy in 1000 Hz - 2000 Hz
    length = (length >> 1);
    WebRtcVad_LogOfEnergy(hp_60, length, kOffsetVector[3], &power,
        &out_vector[3]);

    // Split at 500 Hz
    frequency_band = 3;
    in_ptr = lp_60;
    hp_out_ptr = hp_120;
    lp_out_ptr = lp_120;

    WebRtcVad_SplitFilter(in_ptr, length, &inst->upper_state[frequency_band],
        &inst->lower_state[frequency_band], hp_out_ptr,
        lp_out_ptr);

    // Energy in 500 Hz - 1000 Hz
    length = (length >> 1);
    WebRtcVad_LogOfEnergy(hp_120, length, kOffsetVector[2], &power,
        &out_vector[2]);

    // Split at 250 Hz
    frequency_band = 4;
    in_ptr = lp_120;
    hp_out_ptr = hp_60;
    lp_out_ptr = lp_60;

    WebRtcVad_SplitFilter(in_ptr, length, &inst->upper_state[frequency_band],
        &inst->lower_state[frequency_band], hp_out_ptr,
        lp_out_ptr);

    // Energy in 250 Hz - 500 Hz
    length = (length >> 1);
    WebRtcVad_LogOfEnergy(hp_60, length, kOffsetVector[1], &power,
        &out_vector[1]);

    // Remove DC and LFs
    WebRtcVad_HpOutput(lp_60, length, inst->hp_filter_state, hp_120);

    // Power in 80 Hz - 250 Hz
    WebRtcVad_LogOfEnergy(hp_120, length, kOffsetVector[0], &power,
        &out_vector[0]);

    return power;
}





// Initialize VAD
int WebRtcVad_InitCore(VadInstT *inst, short mode)
{
    int i;

    // Initialization of struct
    inst->vad = 1;
    inst->frame_counter = 0;
    inst->over_hang = 0;
    inst->num_of_speech = 0;

    // Initialization of downsampling filter state
    inst->downsampling_filter_states[0] = 0;
    inst->downsampling_filter_states[1] = 0;
    inst->downsampling_filter_states[2] = 0;
    inst->downsampling_filter_states[3] = 0;

    // Read initial PDF parameters
    for (i = 0; i < NUM_TABLE_VALUES; i++)
    {
        inst->noise_means[i] = kNoiseDataMeans[i];
        inst->speech_means[i] = kSpeechDataMeans[i];
        inst->noise_stds[i] = kNoiseDataStds[i];
        inst->speech_stds[i] = kSpeechDataStds[i];
    }

    // Index and Minimum value vectors are initialized
    for (i = 0; i < 16 * NUM_CHANNELS; i++)
    {
        inst->low_value_vector[i] = 10000;
        inst->index_vector[i] = 0;
    }

    for (i = 0; i < 5; i++)
    {
        inst->upper_state[i] = 0;
        inst->lower_state[i] = 0;
    }

    for (i = 0; i < 4; i++)
    {
        inst->hp_filter_state[i] = 0;
    }

    // Init mean value memory, for FindMin function
    inst->mean_value[0] = 1600;
    inst->mean_value[1] = 1600;
    inst->mean_value[2] = 1600;
    inst->mean_value[3] = 1600;
    inst->mean_value[4] = 1600;
    inst->mean_value[5] = 1600;

    if (mode == 0)
    {
        // Quality mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_Q; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_Q; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_Q; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_Q; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_Q; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_Q; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_Q;
        inst->individual[1] = INDIVIDUAL_20MS_Q;
        inst->individual[2] = INDIVIDUAL_30MS_Q;

        inst->total[0] = TOTAL_10MS_Q;
        inst->total[1] = TOTAL_20MS_Q;
        inst->total[2] = TOTAL_30MS_Q;
    } else if (mode == 1)
    {
        // Low bitrate mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_LBR; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_LBR; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_LBR; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_LBR; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_LBR; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_LBR; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_LBR;
        inst->individual[1] = INDIVIDUAL_20MS_LBR;
        inst->individual[2] = INDIVIDUAL_30MS_LBR;

        inst->total[0] = TOTAL_10MS_LBR;
        inst->total[1] = TOTAL_20MS_LBR;
        inst->total[2] = TOTAL_30MS_LBR;
    } else if (mode == 2)
    {
        // Aggressive mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_AGG; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_AGG; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_AGG; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_AGG; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_AGG; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_AGG; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_AGG;
        inst->individual[1] = INDIVIDUAL_20MS_AGG;
        inst->individual[2] = INDIVIDUAL_30MS_AGG;

        inst->total[0] = TOTAL_10MS_AGG;
        inst->total[1] = TOTAL_20MS_AGG;
        inst->total[2] = TOTAL_30MS_AGG;
    } else
    {
        // Very aggressive mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_VAG; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_VAG; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_VAG; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_VAG; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_VAG; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_VAG; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_VAG;
        inst->individual[1] = INDIVIDUAL_20MS_VAG;
        inst->individual[2] = INDIVIDUAL_30MS_VAG;

        inst->total[0] = TOTAL_10MS_VAG;
        inst->total[1] = TOTAL_20MS_VAG;
        inst->total[2] = TOTAL_30MS_VAG;
    }

    inst->init_flag = kInitCheck;

    return 0;
}

// Set aggressiveness mode
int WebRtcVad_set_mode_core(VadInstT *inst, short mode)
{

    if (mode == 0)
    {
        // Quality mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_Q; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_Q; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_Q; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_Q; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_Q; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_Q; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_Q;
        inst->individual[1] = INDIVIDUAL_20MS_Q;
        inst->individual[2] = INDIVIDUAL_30MS_Q;

        inst->total[0] = TOTAL_10MS_Q;
        inst->total[1] = TOTAL_20MS_Q;
        inst->total[2] = TOTAL_30MS_Q;
    } else if (mode == 1)
    {
        // Low bitrate mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_LBR; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_LBR; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_LBR; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_LBR; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_LBR; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_LBR; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_LBR;
        inst->individual[1] = INDIVIDUAL_20MS_LBR;
        inst->individual[2] = INDIVIDUAL_30MS_LBR;

        inst->total[0] = TOTAL_10MS_LBR;
        inst->total[1] = TOTAL_20MS_LBR;
        inst->total[2] = TOTAL_30MS_LBR;
    } else if (mode == 2)
    {
        // Aggressive mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_AGG; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_AGG; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_AGG; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_AGG; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_AGG; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_AGG; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_AGG;
        inst->individual[1] = INDIVIDUAL_20MS_AGG;
        inst->individual[2] = INDIVIDUAL_30MS_AGG;

        inst->total[0] = TOTAL_10MS_AGG;
        inst->total[1] = TOTAL_20MS_AGG;
        inst->total[2] = TOTAL_30MS_AGG;
    } else if (mode == 3)
    {
        // Very aggressive mode
        inst->over_hang_max_1[0] = OHMAX1_10MS_VAG; // Overhang short speech burst
        inst->over_hang_max_1[1] = OHMAX1_20MS_VAG; // Overhang short speech burst
        inst->over_hang_max_1[2] = OHMAX1_30MS_VAG; // Overhang short speech burst
        inst->over_hang_max_2[0] = OHMAX2_10MS_VAG; // Overhang long speech burst
        inst->over_hang_max_2[1] = OHMAX2_20MS_VAG; // Overhang long speech burst
        inst->over_hang_max_2[2] = OHMAX2_30MS_VAG; // Overhang long speech burst

        inst->individual[0] = INDIVIDUAL_10MS_VAG;
        inst->individual[1] = INDIVIDUAL_20MS_VAG;
        inst->individual[2] = INDIVIDUAL_30MS_VAG;

        inst->total[0] = TOTAL_10MS_VAG;
        inst->total[1] = TOTAL_20MS_VAG;
        inst->total[2] = TOTAL_30MS_VAG;
    } else
    {
        return -1;
    }

    return 0;
}
 
WebRtc_Word16 WebRtcVad_CalcVad32khz(VadInstT *inst, WebRtc_Word16 *speech_frame,
                                     int frame_length)
{
    WebRtc_Word16 len, vad;
    WebRtc_Word16 speechWB[480]; 
    // Downsampled speech frame: 960 samples (30ms in SWB)
    WebRtc_Word16 speechNB[240]; 
    // Downsampled speech frame: 480 samples (30ms in WB)


    // Downsample signal 32->16->8 before doing VAD
    WebRtcVad_Downsampling(speech_frame, speechWB, &(inst->downsampling_filter_states[2]),
                           frame_length);
    len = WEBRTC_SPL_RSHIFT_W16(frame_length, 1);

    WebRtcVad_Downsampling(speechWB, speechNB, inst->downsampling_filter_states, len);
    len = WEBRTC_SPL_RSHIFT_W16(len, 1);

    // Do VAD on an 8 kHz signal
    vad = WebRtcVad_CalcVad8khz(inst, speechNB, len);

    return vad;
}

WebRtc_Word16 WebRtcVad_CalcVad16khz(VadInstT *inst, WebRtc_Word16 *speech_frame,
                                     int frame_length)
{
    WebRtc_Word16 len, vad;
    WebRtc_Word16 speechNB[240]; // Downsampled speech frame: 480 samples (30ms in WB)

    // Wideband: Downsample signal before doing VAD
    WebRtcVad_Downsampling(speech_frame, speechNB, inst->downsampling_filter_states,
                           frame_length);

    len = WEBRTC_SPL_RSHIFT_W16(frame_length, 1);
    vad = WebRtcVad_CalcVad8khz(inst, speechNB, len);

    return vad;
}

WebRtc_Word16 WebRtcVad_CalcVad8khz(VadInstT *inst, WebRtc_Word16 *speech_frame,
                                    int frame_length)
{
    WebRtc_Word16 feature_vector[NUM_CHANNELS], total_power;

    // Get power in the bands
    total_power = WebRtcVad_get_features(inst, speech_frame, frame_length, feature_vector);

    // Make a VAD
    inst->vad = WebRtcVad_GmmProbability(inst, feature_vector, total_power, frame_length);

    return inst->vad;
}

// Calculate probability for both speech and background noise, and perform a
// hypothesis-test.
WebRtc_Word16 WebRtcVad_GmmProbability(VadInstT *inst, WebRtc_Word16 *feature_vector,
                                       WebRtc_Word16 total_power, int frame_length)
{
    int n, k;
    WebRtc_Word16 backval;
    WebRtc_Word16 h0, h1;
    WebRtc_Word16 ratvec, xval;
    WebRtc_Word16 vadflag;
    WebRtc_Word16 shifts0, shifts1;
    WebRtc_Word16 tmp16, tmp16_1, tmp16_2;
    WebRtc_Word16 diff, nr, pos;
    WebRtc_Word16 nmk, nmk2, nmk3, smk, smk2, nsk, ssk;
    WebRtc_Word16 delt, ndelt;
    WebRtc_Word16 maxspe, maxmu;
    WebRtc_Word16 deltaN[NUM_TABLE_VALUES], deltaS[NUM_TABLE_VALUES];
    WebRtc_Word16 ngprvec[NUM_TABLE_VALUES], sgprvec[NUM_TABLE_VALUES];
    WebRtc_Word32 h0test, h1test;
    WebRtc_Word32 tmp32_1, tmp32_2;
    WebRtc_Word32 dotVal;
    WebRtc_Word32 nmid, smid;
    WebRtc_Word32 probn[NUM_MODELS], probs[NUM_MODELS];
    WebRtc_Word16 *nmean1ptr, *nmean2ptr, *smean1ptr, *smean2ptr, *nstd1ptr, *nstd2ptr,
            *sstd1ptr, *sstd2ptr;
    WebRtc_Word16 overhead1, overhead2, individualTest, totalTest;

    // Set the thresholds to different values based on frame length
    if (frame_length == 80)
    {
        // 80 input samples
        overhead1 = inst->over_hang_max_1[0];
        overhead2 = inst->over_hang_max_2[0];
        individualTest = inst->individual[0];
        totalTest = inst->total[0];
    } else if (frame_length == 160)
    {
        // 160 input samples
        overhead1 = inst->over_hang_max_1[1];
        overhead2 = inst->over_hang_max_2[1];
        individualTest = inst->individual[1];
        totalTest = inst->total[1];
    } else
    {
        // 240 input samples
        overhead1 = inst->over_hang_max_1[2];
        overhead2 = inst->over_hang_max_2[2];
        individualTest = inst->individual[2];
        totalTest = inst->total[2];
    }

    if (total_power > MIN_ENERGY)
    { // If signal present at all

        // Set pointers to the gaussian parameters
        nmean1ptr = &inst->noise_means[0];
        nmean2ptr = &inst->noise_means[NUM_CHANNELS];
        smean1ptr = &inst->speech_means[0];
        smean2ptr = &inst->speech_means[NUM_CHANNELS];
        nstd1ptr = &inst->noise_stds[0];
        nstd2ptr = &inst->noise_stds[NUM_CHANNELS];
        sstd1ptr = &inst->speech_stds[0];
        sstd2ptr = &inst->speech_stds[NUM_CHANNELS];

        vadflag = 0;
        dotVal = 0;
        for (n = 0; n < NUM_CHANNELS; n++)
        { // For all channels

            pos = WEBRTC_SPL_LSHIFT_W16(n, 1);
            xval = feature_vector[n];
            //这里，一帧6个子代E；；给了2个模型即H0,H1,分别对应一个二维高斯模型，每个高斯模型有自己均值和方差。
            // Probability for Noise, Q7 * Q20 = Q27
            tmp32_1 = WebRtcVad_GaussianProbability(xval, *nmean1ptr++, *nstd1ptr++,
                                                    &deltaN[pos]);
            probn[0] = (WebRtc_Word32)(kNoiseDataWeights[n] * tmp32_1);
            tmp32_1 = WebRtcVad_GaussianProbability(xval, *nmean2ptr++, *nstd2ptr++,
                                                    &deltaN[pos + 1]);
            probn[1] = (WebRtc_Word32)(kNoiseDataWeights[n + NUM_CHANNELS] * tmp32_1);
            h0test = probn[0] + probn[1]; // Q27
            h0 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(h0test, 12); // Q15

            // Probability for Speech
            tmp32_1 = WebRtcVad_GaussianProbability(xval, *smean1ptr++, *sstd1ptr++,
                                                    &deltaS[pos]);
            probs[0] = (WebRtc_Word32)(kSpeechDataWeights[n] * tmp32_1);
            tmp32_1 = WebRtcVad_GaussianProbability(xval, *smean2ptr++, *sstd2ptr++,
                                                    &deltaS[pos + 1]);
            probs[1] = (WebRtc_Word32)(kSpeechDataWeights[n + NUM_CHANNELS] * tmp32_1);
            h1test = probs[0] + probs[1]; // Q27
            h1 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(h1test, 12); // Q15

            // Get likelihood ratio. Approximate log2(H1/H0) with shifts0 - shifts1
            shifts0 = WebRtcSpl_NormW32(h0test);
            shifts1 = WebRtcSpl_NormW32(h1test);

            if ((h0test > 0) && (h1test > 0))
            {
                ratvec = shifts0 - shifts1;
            } else if (h1test > 0)
            {
                ratvec = 31 - shifts1;
            } else if (h0test > 0)
            {
                ratvec = shifts0 - 31;
            } else
            {
                ratvec = 0;
            }

            // VAD decision with spectrum weighting
            dotVal += WEBRTC_SPL_MUL_16_16(ratvec, kSpectrumWeight[n]);

            // Individual channel test
            if ((ratvec << 2) > individualTest)
            {
                vadflag = 1;
            }

            // Probabilities used when updating model
            if (h0 > 0)
            {
                tmp32_1 = probn[0] & 0xFFFFF000; // Q27
                tmp32_2 = WEBRTC_SPL_LSHIFT_W32(tmp32_1, 2); // Q29
                ngprvec[pos] = (WebRtc_Word16)WebRtcSpl_DivW32W16(tmp32_2, h0);
                ngprvec[pos + 1] = 16384 - ngprvec[pos];
            } else
            {
                ngprvec[pos] = 16384;
                ngprvec[pos + 1] = 0;
            }

            // Probabilities used when updating model
            if (h1 > 0)
            {
                tmp32_1 = probs[0] & 0xFFFFF000;
                tmp32_2 = WEBRTC_SPL_LSHIFT_W32(tmp32_1, 2);
                sgprvec[pos] = (WebRtc_Word16)WebRtcSpl_DivW32W16(tmp32_2, h1);
                sgprvec[pos + 1] = 16384 - sgprvec[pos];
            } else
            {
                sgprvec[pos] = 0;
                sgprvec[pos + 1] = 0;
            }
        }

        // Overall test
        if (dotVal >= totalTest)
        {
            vadflag |= 1;
        }

        // Set pointers to the means and standard deviations.
        nmean1ptr = &inst->noise_means[0];
        smean1ptr = &inst->speech_means[0];
        nstd1ptr = &inst->noise_stds[0];
        sstd1ptr = &inst->speech_stds[0];

        maxspe = 12800;

        // Update the model's parameters
        for (n = 0; n < NUM_CHANNELS; n++)
        {

            pos = WEBRTC_SPL_LSHIFT_W16(n, 1);

            // Get min value in past which is used for long term correction
            backval = WebRtcVad_FindMinimum(inst, feature_vector[n], n); // Q4

            // Compute the "global" mean, that is the sum of the two means weighted
            nmid = WEBRTC_SPL_MUL_16_16(kNoiseDataWeights[n], *nmean1ptr); // Q7 * Q7
            nmid += WEBRTC_SPL_MUL_16_16(kNoiseDataWeights[n+NUM_CHANNELS],
                    *(nmean1ptr+NUM_CHANNELS));
            tmp16_1 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(nmid, 6); // Q8

            for (k = 0; k < NUM_MODELS; k++)
            {

                nr = pos + k;

                nmean2ptr = nmean1ptr + k * NUM_CHANNELS;
                smean2ptr = smean1ptr + k * NUM_CHANNELS;
                nstd2ptr = nstd1ptr + k * NUM_CHANNELS;
                sstd2ptr = sstd1ptr + k * NUM_CHANNELS;
                nmk = *nmean2ptr;
                smk = *smean2ptr;
                nsk = *nstd2ptr;
                ssk = *sstd2ptr;

                // Update noise mean vector if the frame consists of noise only
                nmk2 = nmk;
                if (!vadflag)
                {
                    // deltaN = (x-mu)/sigma^2
                    // ngprvec[k] = probn[k]/(probn[0] + probn[1])

                    delt = (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(ngprvec[nr],
                            deltaN[nr], 11); // Q14*Q11
                    nmk2 = nmk + (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(delt,
                            kNoiseUpdateConst,
                            22); // Q7+(Q14*Q15>>22)
                }

                // Long term correction of the noise mean
                ndelt = WEBRTC_SPL_LSHIFT_W16(backval, 4);
                ndelt -= tmp16_1; // Q8 - Q8
                nmk3 = nmk2 + (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(ndelt,
                        kBackEta,
                        9); // Q7+(Q8*Q8)>>9

                // Control that the noise mean does not drift to much
                tmp16 = WEBRTC_SPL_LSHIFT_W16(k+5, 7);
                if (nmk3 < tmp16)
                    nmk3 = tmp16;
                tmp16 = WEBRTC_SPL_LSHIFT_W16(72+k-n, 7);
                if (nmk3 > tmp16)
                    nmk3 = tmp16;
                *nmean2ptr = nmk3;

                if (vadflag)
                {
                    // Update speech mean vector:
                    // deltaS = (x-mu)/sigma^2
                    // sgprvec[k] = probn[k]/(probn[0] + probn[1])

                    delt = (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(sgprvec[nr],
                            deltaS[nr],
                            11); // (Q14*Q11)>>11=Q14
                    tmp16 = (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(delt,
                            kSpeechUpdateConst,
                            21) + 1;
                    smk2 = smk + (tmp16 >> 1); // Q7 + (Q14 * Q15 >> 22)

                    // Control that the speech mean does not drift to much
                    maxmu = maxspe + 640;
                    if (smk2 < kMinimumMean[k])
                        smk2 = kMinimumMean[k];
                    if (smk2 > maxmu)
                        smk2 = maxmu;

                    *smean2ptr = smk2;

                    // (Q7>>3) = Q4
                    tmp16 = WEBRTC_SPL_RSHIFT_W16((smk + 4), 3);

                    tmp16 = feature_vector[n] - tmp16; // Q4
                    tmp32_1 = WEBRTC_SPL_MUL_16_16_RSFT(deltaS[nr], tmp16, 3);
                    tmp32_2 = tmp32_1 - (WebRtc_Word32)4096; // Q12
                    tmp16 = WEBRTC_SPL_RSHIFT_W16((sgprvec[nr]), 2);
                    tmp32_1 = (WebRtc_Word32)(tmp16 * tmp32_2);// (Q15>>3)*(Q14>>2)=Q12*Q12=Q24

                    tmp32_2 = WEBRTC_SPL_RSHIFT_W32(tmp32_1, 4); // Q20

                    // 0.1 * Q20 / Q7 = Q13
                    if (tmp32_2 > 0)
                        tmp16 = (WebRtc_Word16)WebRtcSpl_DivW32W16(tmp32_2, ssk * 10);
                    else
                    {
                        tmp16 = (WebRtc_Word16)WebRtcSpl_DivW32W16(-tmp32_2, ssk * 10);
                        tmp16 = -tmp16;
                    }
                    // divide by 4 giving an update factor of 0.025
                    tmp16 += 128; // Rounding
                    ssk += WEBRTC_SPL_RSHIFT_W16(tmp16, 8);
                    // Division with 8 plus Q7
                    if (ssk < MIN_STD)
                        ssk = MIN_STD;
                    *sstd2ptr = ssk;
                } else
                {
                    // Update GMM variance vectors
                    // deltaN * (feature_vector[n] - nmk) - 1, Q11 * Q4
                    tmp16 = feature_vector[n] - WEBRTC_SPL_RSHIFT_W16(nmk, 3);

                    // (Q15>>3) * (Q14>>2) = Q12 * Q12 = Q24
                    tmp32_1 = WEBRTC_SPL_MUL_16_16_RSFT(deltaN[nr], tmp16, 3) - 4096;
                    tmp16 = WEBRTC_SPL_RSHIFT_W16((ngprvec[nr]+2), 2);
                    tmp32_2 = (WebRtc_Word32)(tmp16 * tmp32_1);
                    tmp32_1 = WEBRTC_SPL_RSHIFT_W32(tmp32_2, 14);
                    // Q20  * approx 0.001 (2^-10=0.0009766)

                    // Q20 / Q7 = Q13
                    tmp16 = (WebRtc_Word16)WebRtcSpl_DivW32W16(tmp32_1, nsk);
                    if (tmp32_1 > 0)
                        tmp16 = (WebRtc_Word16)WebRtcSpl_DivW32W16(tmp32_1, nsk);
                    else
                    {
                        tmp16 = (WebRtc_Word16)WebRtcSpl_DivW32W16(-tmp32_1, nsk);
                        tmp16 = -tmp16;
                    }
                    tmp16 += 32; // Rounding
                    nsk += WEBRTC_SPL_RSHIFT_W16(tmp16, 6);

                    if (nsk < MIN_STD)
                        nsk = MIN_STD;

                    *nstd2ptr = nsk;
                }
            }

            // Separate models if they are too close - nmid in Q14
            nmid = WEBRTC_SPL_MUL_16_16(kNoiseDataWeights[n], *nmean1ptr);
            nmid += WEBRTC_SPL_MUL_16_16(kNoiseDataWeights[n+NUM_CHANNELS], *nmean2ptr);

            // smid in Q14
            smid = WEBRTC_SPL_MUL_16_16(kSpeechDataWeights[n], *smean1ptr);
            smid += WEBRTC_SPL_MUL_16_16(kSpeechDataWeights[n+NUM_CHANNELS], *smean2ptr);

            // diff = "global" speech mean - "global" noise mean
            diff = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(smid, 9);
            tmp16 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(nmid, 9);
            diff -= tmp16;

            if (diff < kMinimumDifference[n])
            {

                tmp16 = kMinimumDifference[n] - diff; // Q5

                // tmp16_1 = ~0.8 * (kMinimumDifference - diff) in Q7
                // tmp16_2 = ~0.2 * (kMinimumDifference - diff) in Q7
                tmp16_1 = (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(13, tmp16, 2);
                tmp16_2 = (WebRtc_Word16)WEBRTC_SPL_MUL_16_16_RSFT(3, tmp16, 2);

                // First Gauss, speech model
                tmp16 = tmp16_1 + *smean1ptr;
                *smean1ptr = tmp16;
                smid = WEBRTC_SPL_MUL_16_16(tmp16, kSpeechDataWeights[n]);

                // Second Gauss, speech model
                tmp16 = tmp16_1 + *smean2ptr;
                *smean2ptr = tmp16;
                smid += WEBRTC_SPL_MUL_16_16(tmp16, kSpeechDataWeights[n+NUM_CHANNELS]);

                // First Gauss, noise model
                tmp16 = *nmean1ptr - tmp16_2;
                *nmean1ptr = tmp16;

                nmid = WEBRTC_SPL_MUL_16_16(tmp16, kNoiseDataWeights[n]);

                // Second Gauss, noise model
                tmp16 = *nmean2ptr - tmp16_2;
                *nmean2ptr = tmp16;
                nmid += WEBRTC_SPL_MUL_16_16(tmp16, kNoiseDataWeights[n+NUM_CHANNELS]);
            }

            // Control that the speech & noise means do not drift to much
            maxspe = kMaximumSpeech[n];
            tmp16_2 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(smid, 7);
            if (tmp16_2 > maxspe)
            { // Upper limit of speech model
                tmp16_2 -= maxspe;

                *smean1ptr -= tmp16_2;
                *smean2ptr -= tmp16_2;
            }

            tmp16_2 = (WebRtc_Word16)WEBRTC_SPL_RSHIFT_W32(nmid, 7);
            if (tmp16_2 > kMaximumNoise[n])
            {
                tmp16_2 -= kMaximumNoise[n];

                *nmean1ptr -= tmp16_2;
                *nmean2ptr -= tmp16_2;
            }

            nmean1ptr++;
            smean1ptr++;
            nstd1ptr++;
            sstd1ptr++;
        }
        inst->frame_counter++;
    } else
    {
        vadflag = 0;
    }

    // Hangover smoothing
    if (!vadflag)
    {
        if (inst->over_hang > 0)
        {
            vadflag = 2 + inst->over_hang;
            inst->over_hang = inst->over_hang - 1;
        }
        inst->num_of_speech = 0;
    } else
    {
        inst->num_of_speech = inst->num_of_speech + 1;
        if (inst->num_of_speech > NSP_MAX)
        {
            inst->num_of_speech = NSP_MAX;
            inst->over_hang = overhead2;
        } else
            inst->over_hang = overhead1;
    }
    return vadflag;
}

//webrtc_vad.c
