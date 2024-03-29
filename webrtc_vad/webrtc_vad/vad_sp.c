/*
 *  Copyright (c) 2011 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vad_sp.h"

#include <assert.h>

//#include "signal_processing_library.h"
#include "typedefs.h"
//#include "vad_defines.h"
/*


static const int16_t kAllPassCoefsQ13[2] = { 5243, 1392 };  // Q13
#define WEBRTC_SPL_WORD16_MAX       32767
#define WEBRTC_SPL_WORD16_MIN       -32768

// TODO(bjornv): Move this function to vad_filterbank.c.
// Downsampling filter based on splitting filter and allpass functions.
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

// Inserts |feature_value| into |low_value_vector|, if it is one of the 16
// smallest values the last 100 frames. Then calculates and returns the median
// of the five smallest values.
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
*/