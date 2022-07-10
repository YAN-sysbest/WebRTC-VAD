#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "../include/webrtc_vad.h"
 #include "webrtc_vad.h"
 

typedef struct wav_header {
	uint32_t riff_id;
	uint32_t riff_sz;
	uint32_t riff_fmt;
	uint32_t fmt_id;
	uint32_t fmt_sz;
	uint16_t audio_format;
	uint16_t num_channels;
	uint32_t sample_rate;
	uint32_t byte_rate;
	uint16_t block_align;
	uint16_t bits_per_sample;
	uint32_t data_id;
	uint32_t data_sz;
}WAV_HEADER;




#define MAX_FRAMES_CNT_ONESPEECH  100*30 //max length of one speech : 30s
#define MAX_BREAK_MID_SPEECH   50 //1s
#define MIN_SPEECH_LEN         50  //500ms
#define FRAME_NUM 192


typedef struct speechBuffer {
	short buffer[MAX_FRAMES_CNT_ONESPEECH][160];
	short currFrameIndex;
	short speech_len;
	short ns_tail_len;
}SpeechBuffer;

int QuerySpeechBufferIsEmpty(SpeechBuffer* stBuffer) {

	if (stBuffer && stBuffer->currFrameIndex > 0) {
		return 0;
	}

	return 1;

}

int DataInSpeechBuffer(SpeechBuffer* stBuffer, short* Inframe, int is_speech) {

	if (stBuffer->currFrameIndex >= MAX_FRAMES_CNT_ONESPEECH) {
		return -2;//SpeechBuffer is already full
	}

	if (stBuffer && Inframe) {

		memcpy(stBuffer->buffer[stBuffer->currFrameIndex], Inframe, 160 * sizeof(short));
		//2 byte* 160 frames once
		stBuffer->currFrameIndex++;

		if (1 == is_speech) {
			stBuffer->speech_len = stBuffer->currFrameIndex;
			stBuffer->ns_tail_len = 0;

		}
		else {
			stBuffer->ns_tail_len++;

		}
		return 0;
	}
	else {
		return -1;
	}


}



int ClearSpeechBuffer(SpeechBuffer* stBuffer) {
	if (stBuffer) {
		memset(stBuffer->buffer, 0, sizeof(short) * MAX_FRAMES_CNT_ONESPEECH * 160);
		stBuffer->currFrameIndex = 0;
		stBuffer->speech_len = 0;
		stBuffer->ns_tail_len = 0;
	}
	return 0;
}

void TestVAD_my(char* pAudioFile, int nSample, int nMode,char* outAudioFile )
{
	VadInst* pVad = NULL;
	if (WebRtcVad_Create(&pVad))
	{
		perror("WebRtcVad_Create failed!");
		return;
	}

	if (WebRtcVad_Init(pVad))
	{
		perror("WebRtcVad_Init failed!");
		return;
	}

	if (WebRtcVad_set_mode(pVad, nMode))
	{
		perror("WebRtcVad_set_mode failed!");
		return;
	}

	FILE* fp = NULL;
	FILE* fw ;
	FILE* fs;
	FILE* fpR = NULL;
	int i;
	fopen_s(&fp, pAudioFile, "rb");
	if (!fp) {
		printf("Error:file \"%s\" not found. ", pAudioFile);
		return;
	}

	fopen_s(&fw, outAudioFile, "wb");
	if (!fw) {
		printf("Error:file \"%s\" not found. ", outAudioFile);
		return;
	}
	
	//fseek(fp, 0, SEEK_END);
	//unsigned int nLen = ftell(fp);
	//fseek(fp, 0, SEEK_SET);
	short shBufferIn[FRAME_NUM] = { 0 };
	short shBuffer400[FRAME_NUM] = { 0 };
	short shBuffer0[FRAME_NUM] = { 0 };

	for (i = 0; i < FRAME_NUM; i++) { shBuffer400[i] = 400; shBuffer0[i] = 0; };


	//short buf_wav_header[100] = { 0 };
	//short len_header = fread(buf_wav_header, 1, sizeof(WAV_HEADER), fp);
	//printf("len_header = %d\n", len_header);


	SpeechBuffer my_buffer;
	ClearSpeechBuffer(&my_buffer);
	int speech_num = 0;
	int frame = 0;
	char ret_file_name[50] = { 0 };

	fopen_s(&fs, "test--realwav.txt", "w");
	if (fs == NULL) { printf("Cann't open the file: %s !!!\n", "test.txt");  exit(0); }

	while (1)
	{
		if (FRAME_NUM != fread(shBufferIn, 2, FRAME_NUM, fp)) break;
		if (frame  % 100 == 0) printf("%ld points have been processed...\n", frame * FRAME_NUM);
		frame++;
		int nRet = WebRtcVad_Process(pVad, 16000, shBufferIn, FRAME_NUM);
		fprintf_s(fs, "%d\n", nRet);
		if (nRet == 1)//检测到为人声
		{
			fwrite(shBufferIn, 2, FRAME_NUM, fw);
			//fwrite(shBuffer0, 2, 160, fw);

		} 
		else { //检测到非人声
			fwrite(shBuffer0, 2, FRAME_NUM, fw);
			//fwrite(shBufferIn, 2, 160, fw);
		} 
	}

	fclose(fp); fclose(fs); fclose(fw);
	WebRtcVad_Free(pVad);
}



int main()
{ 
	TestVAD_my("sp0102noise.pcm", 16000, 3, "out.pcm");
	return 0;
}

