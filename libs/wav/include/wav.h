#ifndef WAV_INCLUDED
#define WAV_INCLUDED

#include <stdint.h>

/* PCM type - Regardless of the bit depth of the file, all data is handled in memory as signed 32-bit data. */
typedef int32_t WAVPcmData;

/* WAV data format */
typedef enum WAVDataFormatTag {
    WAV_DATA_FORMAT_PCM             /* PCM only */
} WAVDataFormat;

/* API result type */
typedef enum WAVApiResultTag {
    WAV_APIRESULT_OK = 0,
    WAV_APIRESULT_NG,
    WAV_APIRESULT_INVALID_FORMAT,     /* Invalid format */
    WAV_APIRESULT_IOERROR,            /* File I/O error */
    WAV_APIRESULT_INVALID_PARAMETER   /* invalid argument */
} WAVApiResult;

/* WAV file format */
struct WAVFileFormat {
    WAVDataFormat data_format;      /* Data format */
    uint32_t      num_channels;     /* Number of channels */
    uint32_t      sampling_rate;    /* Sampling rate */
    uint32_t      bits_per_sample;  /* Quantization bit rate */
    uint32_t      num_samples;      /* Number of samples */
};

/* WAV file handle */
struct WAVFile {
    struct WAVFileFormat  format;   /* Format */
    WAVPcmData**          data;     /* Actual data */
};

/* accessor */
#define WAVFile_PCM(wavfile, samp, ch)  (wavfile->data[(ch)][(samp)])

#ifdef __cplusplus
extern "C" {
#endif

/* Create a WAV file handle from the file */
struct WAVFile* WAV_CreateFromFile(const char* filename);

/* Create a new WAV file handle with the specified format */
struct WAVFile* WAV_Create(const struct WAVFileFormat* format);

/* Discard the WAV file handle */
void WAV_Destroy(struct WAVFile* wavfile);

/* Write file */
WAVApiResult WAV_WriteToFile(
        const char* filename, const struct WAVFile* wavfile);

/* Read only WAV file formats from file */
WAVApiResult WAV_GetWAVFormatFromFile(
        const char* filename, struct WAVFileFormat* format);

#ifdef __cplusplus
}
#endif

#endif /* WAV_INCLUDED */
