 
#ifndef WEBRTC_VAD_WEBRTC_VAD_H_
#define WEBRTC_VAD_WEBRTC_VAD_H_

#include "typedefs.h"

typedef struct WebRtcVadInst VadInst;

#ifdef __cplusplus
extern "C"
{
#endif 


 
WebRtc_Word16 WebRtcVad_Create(VadInst **vad_inst);
 
WebRtc_Word16 WebRtcVad_Free(VadInst *vad_inst);
 
WebRtc_Word16 WebRtcVad_Init(VadInst *vad_inst);

 
WebRtc_Word16 WebRtcVad_set_mode(VadInst *vad_inst, WebRtc_Word16 mode);
 
WebRtc_Word16 WebRtcVad_Process(VadInst *vad_inst,
                                WebRtc_Word16 fs,
                                WebRtc_Word16 *speech_frame,
                                WebRtc_Word16 frame_length);

#ifdef __cplusplus
}
#endif

#endif // WEBRTC_VAD_WEBRTC_VAD_H_
