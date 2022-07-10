 
#ifndef WEBRTC_TYPEDEFS_H_
#define WEBRTC_TYPEDEFS_H_
 
//define WEBRTC_EXTERN extern
#define G_CONST const
//define WEBRTC_INLINE extern __inline
  
    
    typedef signed char         int8_t;
    typedef signed short        int16_t;
    typedef signed int          int32_t;
    typedef signed long long    int64_t;
    typedef unsigned char       uint8_t;
    typedef unsigned short      uint16_t;
    typedef unsigned int        uint32_t;
    typedef unsigned long long  uint64_t;
 
 
    typedef int32_t             WebRtc_Word32;
    typedef uint32_t            WebRtc_UWord32;
    typedef int16_t             WebRtc_Word16;
    typedef uint16_t            WebRtc_UWord16;
    typedef char                WebRtc_Word8;
    typedef uint8_t             WebRtc_UWord8;

    // Define endian for the platform
    #define WEBRTC_LITTLE_ENDIAN

    //vad_defines

#define NUM_CHANNELS        6   // Eight frequency bands
#define NUM_MODELS          2   // Number of Gaussian models
#define NUM_TABLE_VALUES    NUM_CHANNELS * NUM_MODELS

#define MIN_ENERGY          10
#define ALPHA1              6553    // 0.2 in Q15
#define ALPHA2              32439   // 0.99 in Q15
#define NSP_MAX             6       // Maximum number of VAD=1 frames in a row counted
#define MIN_STD             384     // Minimum standard deviation
// Mode 0, Quality thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_Q   24
#define INDIVIDUAL_20MS_Q   21      // (log10(2)*66)<<2 ~=16
#define INDIVIDUAL_30MS_Q   24

#define TOTAL_10MS_Q        57
#define TOTAL_20MS_Q        48
#define TOTAL_30MS_Q        57

#define OHMAX1_10MS_Q       8  // Max Overhang 1
#define OHMAX2_10MS_Q       14 // Max Overhang 2
#define OHMAX1_20MS_Q       4  // Max Overhang 1
#define OHMAX2_20MS_Q       7  // Max Overhang 2
#define OHMAX1_30MS_Q       3
#define OHMAX2_30MS_Q       5

// Mode 1, Low bitrate thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_LBR 37
#define INDIVIDUAL_20MS_LBR 32
#define INDIVIDUAL_30MS_LBR 37

#define TOTAL_10MS_LBR      100
#define TOTAL_20MS_LBR      80
#define TOTAL_30MS_LBR      100

#define OHMAX1_10MS_LBR     8  // Max Overhang 1
#define OHMAX2_10MS_LBR     14 // Max Overhang 2
#define OHMAX1_20MS_LBR     4
#define OHMAX2_20MS_LBR     7

#define OHMAX1_30MS_LBR     3
#define OHMAX2_30MS_LBR     5

// Mode 2, Very aggressive thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_AGG 82
#define INDIVIDUAL_20MS_AGG 78
#define INDIVIDUAL_30MS_AGG 82

#define TOTAL_10MS_AGG      285 //580
#define TOTAL_20MS_AGG      260
#define TOTAL_30MS_AGG      285

#define OHMAX1_10MS_AGG     6  // Max Overhang 1
#define OHMAX2_10MS_AGG     9  // Max Overhang 2
#define OHMAX1_20MS_AGG     3
#define OHMAX2_20MS_AGG     5

#define OHMAX1_30MS_AGG     2
#define OHMAX2_30MS_AGG     3

// Mode 3, Super aggressive thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_VAG 94
#define INDIVIDUAL_20MS_VAG 94
#define INDIVIDUAL_30MS_VAG 94

#define TOTAL_10MS_VAG      1100 //1700
#define TOTAL_20MS_VAG      1050
#define TOTAL_30MS_VAG      1100

#define OHMAX1_10MS_VAG     6  // Max Overhang 1
#define OHMAX2_10MS_VAG     9  // Max Overhang 2
#define OHMAX1_20MS_VAG     3
#define OHMAX2_20MS_VAG     5

#define OHMAX1_30MS_VAG     2
#define OHMAX2_30MS_VAG     3
     
#endif
     
    //inline 
    /*
    static   int16_t WebRtcSpl_GetSizeInBits(uint32_t n) {
        int bits;

        if (0xFFFF0000 & n) {
            bits = 16;
        }
        else {
            bits = 0;
        }
        if (0x0000FF00 & (n >> bits)) bits += 8;
        if (0x000000F0 & (n >> bits)) bits += 4;
        if (0x0000000C & (n >> bits)) bits += 2;
        if (0x00000002 & (n >> bits)) bits += 1;
        if (0x00000001 & (n >> bits)) bits += 1;

        return bits;
    }

    static   int WebRtcSpl_NormW32(int32_t a) {
        int zeros;

        if (a <= 0) a ^= 0xFFFFFFFF;

        if (!(0xFFFF8000 & a)) {
            zeros = 16;
        }
        else {
            zeros = 0;
        }
        if (!(0xFF800000 & (a << zeros))) zeros += 8;
        if (!(0xF8000000 & (a << zeros))) zeros += 4;
        if (!(0xE0000000 & (a << zeros))) zeros += 2;
        if (!(0xC0000000 & (a << zeros))) zeros += 1;

        return zeros;
    }


    static   int WebRtcSpl_NormU32(uint32_t a) {
        int zeros;

        if (a == 0) return 0;

        if (!(0xFFFF0000 & a)) {
            zeros = 16;
        }
        else {
            zeros = 0;
        }
        if (!(0xFF000000 & (a << zeros))) zeros += 8;
        if (!(0xF0000000 & (a << zeros))) zeros += 4;
        if (!(0xC0000000 & (a << zeros))) zeros += 2;
        if (!(0x80000000 & (a << zeros))) zeros += 1;

        return zeros;
    }


    */




    /*
    
#ifndef WEBRTC_TYPEDEFS_H_
#define WEBRTC_TYPEDEFS_H_
 
#define WEBRTC_EXTERN extern
#define G_CONST const
#define WEBRTC_INLINE extern __inline
 
#if defined(WIN32) 
#elif defined(__APPLE__) 
#else
    // Linux etc.
    #if !defined(WEBRTC_TARGET_PC)
        #define WEBRTC_TARGET_PC
    #endif
#endif
 
#if defined(_M_X64) || defined(__x86_64__)
#define WEBRTC_ARCH_X86_FAMILY
#define WEBRTC_ARCH_X86_64
#define WEBRTC_ARCH_64_BITS
#define WEBRTC_ARCH_LITTLE_ENDIAN
#elif defined(_M_IX86) || defined(__i386__)
#define WEBRTC_ARCH_X86_FAMILY
#define WEBRTC_ARCH_X86
#define WEBRTC_ARCH_32_BITS
#define WEBRTC_ARCH_LITTLE_ENDIAN
#elif defined(__ARMEL__)
// TODO(andrew): We'd prefer to control platform defines here, but this is
// currently provided by the Android makefiles. Commented to avoid duplicate
// definition warnings.
//#define WEBRTC_ARCH_ARM
// TODO(andrew): Chromium uses the following two defines. Should we switch?
//#define WEBRTC_ARCH_ARM_FAMILY
//#define WEBRTC_ARCH_ARMEL
#define WEBRTC_ARCH_32_BITS
#define WEBRTC_ARCH_LITTLE_ENDIAN
#elif defined(__mips__)
#define WEBRTC_ARCH_32_BITS
#define WEBRTC_ARCH_LITTLE_ENDIAN
#else
#error Please add support for your architecture in typedefs.h
#endif

#if defined(__SSE2__) || defined(_MSC_VER)
#define WEBRTC_USE_SSE2
#endif

#if defined(WEBRTC_TARGET_PC)

#if !defined(_MSC_VER)
  #include <stdint.h>
#else 
    
    typedef signed char         int8_t;
    typedef signed short        int16_t;
    typedef signed int          int32_t;
    typedef signed long long    int64_t;
    typedef unsigned char       uint8_t;
    typedef unsigned short      uint16_t;
    typedef unsigned int        uint32_t;
    typedef unsigned long long  uint64_t;
#endif
 
    typedef int32_t             WebRtc_Word32;
    typedef uint32_t            WebRtc_UWord32;
    typedef int16_t             WebRtc_Word16;
    typedef uint16_t            WebRtc_UWord16;
    typedef char                WebRtc_Word8;
    typedef uint8_t             WebRtc_UWord8;

    // Define endian for the platform
    #define WEBRTC_LITTLE_ENDIAN

    //vad_defines

#define NUM_CHANNELS        6   // Eight frequency bands
#define NUM_MODELS          2   // Number of Gaussian models
#define NUM_TABLE_VALUES    NUM_CHANNELS * NUM_MODELS

#define MIN_ENERGY          10
#define ALPHA1              6553    // 0.2 in Q15
#define ALPHA2              32439   // 0.99 in Q15
#define NSP_MAX             6       // Maximum number of VAD=1 frames in a row counted
#define MIN_STD             384     // Minimum standard deviation
// Mode 0, Quality thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_Q   24
#define INDIVIDUAL_20MS_Q   21      // (log10(2)*66)<<2 ~=16
#define INDIVIDUAL_30MS_Q   24

#define TOTAL_10MS_Q        57
#define TOTAL_20MS_Q        48
#define TOTAL_30MS_Q        57

#define OHMAX1_10MS_Q       8  // Max Overhang 1
#define OHMAX2_10MS_Q       14 // Max Overhang 2
#define OHMAX1_20MS_Q       4  // Max Overhang 1
#define OHMAX2_20MS_Q       7  // Max Overhang 2
#define OHMAX1_30MS_Q       3
#define OHMAX2_30MS_Q       5

// Mode 1, Low bitrate thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_LBR 37
#define INDIVIDUAL_20MS_LBR 32
#define INDIVIDUAL_30MS_LBR 37

#define TOTAL_10MS_LBR      100
#define TOTAL_20MS_LBR      80
#define TOTAL_30MS_LBR      100

#define OHMAX1_10MS_LBR     8  // Max Overhang 1
#define OHMAX2_10MS_LBR     14 // Max Overhang 2
#define OHMAX1_20MS_LBR     4
#define OHMAX2_20MS_LBR     7

#define OHMAX1_30MS_LBR     3
#define OHMAX2_30MS_LBR     5

// Mode 2, Very aggressive thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_AGG 82
#define INDIVIDUAL_20MS_AGG 78
#define INDIVIDUAL_30MS_AGG 82

#define TOTAL_10MS_AGG      285 //580
#define TOTAL_20MS_AGG      260
#define TOTAL_30MS_AGG      285

#define OHMAX1_10MS_AGG     6  // Max Overhang 1
#define OHMAX2_10MS_AGG     9  // Max Overhang 2
#define OHMAX1_20MS_AGG     3
#define OHMAX2_20MS_AGG     5

#define OHMAX1_30MS_AGG     2
#define OHMAX2_30MS_AGG     3

// Mode 3, Super aggressive thresholds - Different thresholds for the different frame lengths
#define INDIVIDUAL_10MS_VAG 94
#define INDIVIDUAL_20MS_VAG 94
#define INDIVIDUAL_30MS_VAG 94

#define TOTAL_10MS_VAG      1100 //1700
#define TOTAL_20MS_VAG      1050
#define TOTAL_30MS_VAG      1100

#define OHMAX1_10MS_VAG     6  // Max Overhang 1
#define OHMAX2_10MS_VAG     9  // Max Overhang 2
#define OHMAX1_20MS_VAG     3
#define OHMAX2_20MS_VAG     5

#define OHMAX1_30MS_VAG     2
#define OHMAX2_30MS_VAG     3
     
#endif

#endif  // WEBRTC_TYPEDEFS_H_




    
    
    */