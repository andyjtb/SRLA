#include "wav.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Parser read buffer size */
#define WAVBITBUFFER_BUFFER_SIZE         (128 * 1024)

/* Get the lowest n_bits */
/* Supplementary Note) ((1 << n_bits) - 1) is a mask that extracts only the lower numbers */
#define WAV_GetLowerBits(n_bits, val) ((val) & (uint32_t)((1 << (n_bits)) - 1))

/* Get the smaller value of a or b */
#define WAV_Min(a, b) (((a) < (b)) ? (a) : (b))

/* Internal error type */
typedef enum WAVErrorTag {
    WAV_ERROR_OK = 0,             /* OK */
    WAV_ERROR_NG,                 /* Unclassifiable failure */
    WAV_ERROR_IO,                 /* I/O error */
    WAV_ERROR_INVALID_PARAMETER,  /* Invalid argument */
    WAV_ERROR_INVALID_FORMAT      /* Invalid format */
} WAVError;

/* bit buffer */
struct WAVBitBuffer {
    uint8_t   bytes[WAVBITBUFFER_BUFFER_SIZE];   /* bit buffer */
    uint32_t  bit_count;                        /* bit input count */
    int32_t   byte_pos;                         /* Byte string read position */
};

/* Parser */
struct WAVParser {
    FILE*               fp;       /* read file pointer */
    struct WAVBitBuffer buffer;   /* bit buffer */
};

/* Writer */
struct WAVWriter {
    FILE*     fp;                 /* write file pointer */
    uint32_t  bit_buffer;         /* Bit in the middle of output */
    uint32_t  bit_count;          /* Output count */
    struct WAVBitBuffer buffer;   /* bit buffer */
};

/* Initialize the parser */
static void WAVParser_Initialize(struct WAVParser* parser, FILE* fp);
/* End of parser use */
static void WAVParser_Finalize(struct WAVParser* parser);
/* Get n_bits and right-justify the result */
static WAVError WAVParser_GetBits(struct WAVParser* parser, uint32_t n_bits, uint64_t* bitsbuf);
/* Seek (compliant with fseek) */
static WAVError WAVParser_Seek(struct WAVParser* parser, int32_t offset, int32_t wherefrom);
/* Writer initialization */
static void WAVWriter_Initialize(struct WAVWriter* writer, FILE* fp);
/* Writer termination */
static void WAVWriter_Finalize(struct WAVWriter* writer);
/* Write the lowest n_bits of val */
static WAVError WAVWriter_PutBits(struct WAVWriter* writer, uint64_t val, uint32_t n_bits);
/* Clear the bits accumulated in the buffer */
static WAVError WAVWriter_Flush(struct WAVWriter* writer);
/* Output bit pattern in little endian */
static WAVError WAVWriter_PutLittleEndianBytes(
        struct WAVWriter* writer, uint32_t nbytes, uint64_t data);

/* Use a writer to output the header according to the file format */
static WAVError WAVWriter_PutWAVHeader(
        struct WAVWriter* writer, const struct WAVFileFormat* format);
/* Output PCM data using the writer */
static WAVError WAVWriter_PutWAVPcmData(
        struct WAVWriter* writer, const struct WAVFile* wavfile);

/* Get bit pattern in little endian */
static WAVError WAVParser_GetLittleEndianBytes(
        struct WAVParser* parser, uint32_t nbytes, uint64_t* bitsbuf);
/* Get string using parser */
static WAVError WAVParser_GetString(
        struct WAVParser* parser, char* string_buffer, uint32_t string_length);
/* Use parser to get string/check for match */
static WAVError WAVParser_CheckSignatureString(
        struct WAVParser* parser, const char* signature, uint32_t signature_length);
/* Use the parser to read the file format */
static WAVError WAVParser_GetWAVFormat(
        struct WAVParser* parser, struct WAVFileFormat* format);
/* Read the PCM data using the parser */
static WAVError WAVParser_GetWAVPcmData(
        struct WAVParser* parser, struct WAVFile* wavfile);

/* Convert 8bit PCM format to 32bit format */
static int32_t WAV_Convert8bitPCMto32bitPCM(int32_t in_8bitpcm);
/* Convert 16bit PCM format to 32bit format */
static int32_t WAV_Convert16bitPCMto32bitPCM(int32_t in_16bitpcm);
/* Convert 24bit PCM format to 32bit format */
static int32_t WAV_Convert24bitPCMto32bitPCM(int32_t in_24bitpcm);
/* Convert 32bit PCM format to 32bit format */
static int32_t WAV_Convert32bitPCMto32bitPCM(int32_t in_32bitpcm);

/* Convert 32bit PCM format to 32bit format */
static int32_t WAV_Convert32bitPCMto32bitPCM(int32_t in_32bitpcm);

/* Use the parser to read the file format */
static WAVError WAVParser_GetWAVFormat(
        struct WAVParser* parser, struct WAVFileFormat* format)
{
    uint64_t  bitsbuf;
    int32_t   fmt_chunk_size;
    struct WAVFileFormat tmp_format;

    /* Argument check */
    if (parser == NULL || format == NULL) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /* Check for headers 'R', 'I', 'F', 'F' */
    if (WAVParser_CheckSignatureString(parser, "RIFF", 4) != WAV_ERROR_OK) {
        return WAV_ERROR_INVALID_FORMAT;
    }

    /* File size - 8 (skip) */
    if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Check for headers 'W', 'A', 'V', 'E' */
    if (WAVParser_CheckSignatureString(parser, "WAVE", 4) != WAV_ERROR_OK) {
        return WAV_ERROR_INVALID_FORMAT;
    }

    /* Check for fmt chunk headers 'f', 'm', 't', ' ' */
    if (WAVParser_CheckSignatureString(parser, "fmt ", 4) != WAV_ERROR_OK) {
        return WAV_ERROR_INVALID_FORMAT;
    }

    /*
/* Get the number of bytes in the fmt chunk
* Supplementary note/Note) The contents (extensions) of fmt chunks with a size greater than 16 will be skipped. */
*/
    if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    fmt_chunk_size = (int32_t)bitsbuf;

    /*
/* Check the format ID
* Note) Only 1 (Linear PCM) is supported */
*/
    if (WAVParser_GetLittleEndianBytes(parser, 2, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    if (bitsbuf != 1) {
        /* fprintf(stderr, "Unsupported format: fmt chunk format ID \n"); */
        return WAV_ERROR_INVALID_FORMAT;
    }
    tmp_format.data_format = WAV_DATA_FORMAT_PCM;

    /* Number of channels */
    if (WAVParser_GetLittleEndianBytes(parser, 2, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    tmp_format.num_channels = (uint32_t)bitsbuf;

    /* Sampling rate */
    if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    tmp_format.sampling_rate =(uint32_t) bitsbuf;

    /* Data speed (byte/sec) skipped */
    if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Size number per block is skipped */
    if (WAVParser_GetLittleEndianBytes(parser, 2, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Quantization bit depth (bits per sample) */
    if (WAVParser_GetLittleEndianBytes(parser, 2, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    tmp_format.bits_per_sample = (uint32_t)bitsbuf;

    /* Reading of extensions is not supported: Skip reading */
    if (fmt_chunk_size > 16) {
        fprintf(stderr, "Warning: skip fmt chunk extention (unsupported). \n");
        if (WAVParser_Seek(parser, fmt_chunk_size - 16, SEEK_CUR) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    }

    /* Chunk read */
    while (1) {
        char string_buf[4];
        /* Get chunk string */
        if (WAVParser_GetString(parser, string_buf, 4) != WAV_ERROR_OK) {
            return WAV_ERROR_IO;
        }
        if (strncmp(string_buf, "data", 4) == 0) {
            /* End when data chunk is found */
            break;
        } else {
            /* Other chunks are skipped by seeking after getting only their size */
            if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) {
                return WAV_ERROR_IO;
            }
            /* printf("chunk:%s size:%d \n", string_buf, (int32_t)bitsbuf); */
            WAVParser_Seek(parser, (int32_t)bitsbuf, SEEK_CUR);
        }
    }

    /* Number of samples: Calculated from number of waveform data bytes */
    if (WAVParser_GetLittleEndianBytes(parser, 4, &bitsbuf) != WAV_ERROR_OK) { return WAV_ERROR_IO; }
    tmp_format.num_samples = (uint32_t)bitsbuf;
    assert(tmp_format.num_samples % ((tmp_format.bits_per_sample / 8) * tmp_format.num_channels) == 0);
    tmp_format.num_samples /= ((tmp_format.bits_per_sample / 8) * tmp_format.num_channels);

    /* Copy structure */
    *format = tmp_format;

    return WAV_ERROR_OK;
}

/* Read the PCM data using the parser */
static WAVError WAVParser_GetWAVPcmData(
        struct WAVParser* parser, struct WAVFile* wavfile)
{
    uint32_t  ch, sample, bytes_per_sample;
    uint64_t  bitsbuf;
    int32_t   (*convert_to_sint32_func)(int32_t);

    /* Argument check */
    if (parser == NULL || wavfile == NULL) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /* Determine the PCM data conversion function according to the bit depth */
    switch (wavfile->format.bits_per_sample) {
    case 8:
        convert_to_sint32_func = WAV_Convert8bitPCMto32bitPCM;
        break;
    case 16:
        convert_to_sint32_func = WAV_Convert16bitPCMto32bitPCM;
        break;
    case 24:
        convert_to_sint32_func = WAV_Convert24bitPCMto32bitPCM;
        break;
    case 32:
        convert_to_sint32_func = WAV_Convert32bitPCMto32bitPCM;
        break;
    default:
        /* fprintf(stderr, "Unsupported bits per sample format(=%d). \n", wavfile->format.bits_per_sample); */
        return WAV_ERROR_INVALID_FORMAT;
    }

    /* Read data */
    bytes_per_sample = wavfile->format.bits_per_sample / 8;
    for (sample = 0; sample < wavfile->format.num_samples; sample++) {
        for (ch = 0; ch < wavfile->format.num_channels; ch++) {
            if (WAVParser_GetLittleEndianBytes(parser, bytes_per_sample, &bitsbuf) != WAV_ERROR_OK) {
                return WAV_ERROR_IO;
            }
            /* Convert to 32-bit integer format and set to data */
            wavfile->data[ch][sample] = convert_to_sint32_func((int32_t)(bitsbuf));
        }
    }

    return WAV_ERROR_OK;
}

/* Read only WAV file formats from file */
WAVApiResult WAV_GetWAVFormatFromFile(
        const char* filename, struct WAVFileFormat* format)
{
    struct WAVParser parser;
    FILE*            fp;

    /* Argument check */
    if (filename == NULL || format == NULL) {
        return WAV_APIRESULT_NG;
    }

    /* Open wav file */
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        /* fprintf(stderr, "Failed to open %s. \n", filename); */
        return WAV_APIRESULT_NG;
    }

    /* Parser initialization */
    WAVParser_Initialize(&parser, fp);

    /* Read header */
    if (WAVParser_GetWAVFormat(&parser, format) != WAV_ERROR_OK) {
        return WAV_APIRESULT_NG;
    }

    /* End of parser use */
    WAVParser_Finalize(&parser);

    /* Close the file */
    fclose(fp);

    return WAV_APIRESULT_OK;
}

/* Create a WAV file handle from the file */
struct WAVFile* WAV_CreateFromFile(const char* filename)
{
    struct WAVParser      parser;
    FILE*                 fp;
    struct WAVFile*       wavfile;
    struct WAVFileFormat  format;

    /* Argument check */
    if (filename == NULL) {
        return NULL;
    }

    /* Open wav file */
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        /* fprintf(stderr, "Failed to open %s. \n", filename); */
        return NULL;
    }

    /* Parser initialization */
    WAVParser_Initialize(&parser, fp);

    /* Read header */
    if (WAVParser_GetWAVFormat(&parser, &format) != WAV_ERROR_OK) {
        return NULL;
    }

    /* Create handle */
    wavfile = WAV_Create(&format);
    if (wavfile == NULL) {
        return NULL;
    }

    /* Read PCM data */
    if (WAVParser_GetWAVPcmData(&parser, wavfile) != WAV_ERROR_OK) {
        goto EXIT_FAILURE_WITH_DATA_RELEASE;
    }

    /* Parser end */
    WAVParser_Finalize(&parser);

    /* Close the file */
    fclose(fp);

    /* normal termination */
    return wavfile;

    /* Release all data allocated by the handle and exit */
EXIT_FAILURE_WITH_DATA_RELEASE:
    WAV_Destroy(wavfile);
    WAVParser_Finalize(&parser);
    fclose(fp);
    return NULL;
}

/* Create a new WAV file handle with the specified format */
struct WAVFile* WAV_Create(const struct WAVFileFormat* format)
{
    uint32_t ch;
    struct WAVFile* wavfile;

    /* Argument check */
    if (format == NULL) {
        return NULL;
    }

    /* Currently only PCM format is supported */
    if (format->data_format != WAV_DATA_FORMAT_PCM) {
        /* fprintf(stderr, "Unsupported wav data format. \n"); */
        return NULL;
    }

    /* Create handle */
    wavfile = (struct WAVFile *)malloc(sizeof(struct WAVFile));
    if (wavfile == NULL) {
        goto EXIT_FAILURE_WITH_DATA_RELEASE;
    }

    /* Get format information by copying the structure */
    wavfile->format = (*format);

    /* Data area allocation */
    wavfile->data = (WAVPcmData **)malloc(sizeof(WAVPcmData *) * format->num_channels);
    if (wavfile->data == NULL) {
        goto EXIT_FAILURE_WITH_DATA_RELEASE;
    }
    for (ch = 0; ch < format->num_channels; ch++) {
        wavfile->data[ch] = (WAVPcmData *)calloc(format->num_samples, sizeof(WAVPcmData));
        if (wavfile->data[ch] == NULL) {
            goto EXIT_FAILURE_WITH_DATA_RELEASE;
        }
    }

    return wavfile;

EXIT_FAILURE_WITH_DATA_RELEASE:
    WAV_Destroy(wavfile);
    return NULL;
}

/* Convert 8bit PCM format to 32bit format */
static int32_t WAV_Convert8bitPCMto32bitPCM(int32_t in_8bitpcm)
{
    /* Subtract 128, which corresponds to silence */
    return (in_8bitpcm - 128);
}

/* Convert 16bit PCM format to 32bit format */
static int32_t WAV_Convert16bitPCMto32bitPCM(int32_t in_16bitpcm)
{
    /* First make it 32-bit wide, then use arithmetic right shift to add the sign */
    return (in_16bitpcm << 16) >> 16;
}

/* Convert 24bit PCM format to 32bit format */
static int32_t WAV_Convert24bitPCMto32bitPCM(int32_t in_24bitpcm)
{
    /* First make it 32-bit wide, then use arithmetic right shift to add the sign */
    return (in_24bitpcm << 8) >> 8;
}

/* Convert 32bit PCM format to 32bit format */
static int32_t WAV_Convert32bitPCMto32bitPCM(int32_t in_32bitpcm)
{
    /* Do nothing */
    return in_32bitpcm;
}

/* Initialize the parser */
static void WAVParser_Initialize(struct WAVParser* parser, FILE* fp)
{
    parser->fp                = fp;
    memset(&parser->buffer, 0, sizeof(struct WAVBitBuffer));
    parser->buffer.byte_pos   = -1;
}

/* End of parser use */
static void WAVParser_Finalize(struct WAVParser* parser)
{
    parser->fp                = NULL;
    memset(&parser->buffer, 0, sizeof(struct WAVBitBuffer));
    parser->buffer.byte_pos   = -1;
}

/* Get n_bits and right-justify the result */
static WAVError WAVParser_GetBits(struct WAVParser* parser, uint32_t n_bits, uint64_t* bitsbuf)
{
    uint64_t tmp;
    struct WAVBitBuffer *buf = &(parser->buffer);

    /* Argument check */
    if (parser == NULL || bitsbuf == NULL || n_bits > 64) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /* First load */
    if (buf->byte_pos == -1) {
        if (fread(buf->bytes, sizeof(uint8_t), WAVBITBUFFER_BUFFER_SIZE, parser->fp) == 0) {
            return WAV_ERROR_IO;
        }
        buf->byte_pos   = 0;
        buf->bit_count  = 8;
    }

    /*
/* Fill data from the most significant bit
* In the first loop, set it to the higher bits of tmp
* From the second loop onwards, input in 8-bit units and set it to tmp */
*/
    tmp = 0;
    while (n_bits > buf->bit_count) {
        /* Fill from the most significant bit */
        n_bits  -= buf->bit_count;
        tmp     |= (uint64_t)WAV_GetLowerBits(buf->bit_count, buf->bytes[buf->byte_pos]) << n_bits;

        /* Read forward 1 byte */
        buf->byte_pos++;
        buf->bit_count   = 8;

        /* If the buffer is full, re-read it */
        if (buf->byte_pos == WAVBITBUFFER_BUFFER_SIZE) {
            if (fread(buf->bytes, sizeof(uint8_t), WAVBITBUFFER_BUFFER_SIZE, parser->fp) == 0) {
                return WAV_ERROR_IO;
            }
            buf->byte_pos = 0;
        }
    }

    /*
/* Processing of fractional bits
* Set the remaining bits to the most significant bit of tmp */
*/
    buf->bit_count -= n_bits;
    tmp            |= (uint64_t)WAV_GetLowerBits(n_bits, (uint32_t)(buf->bytes[buf->byte_pos] >> buf->bit_count));

    *bitsbuf = tmp;
    return WAV_ERROR_OK;
}

/* Seek (compliant with fseek) */
static WAVError WAVParser_Seek(struct WAVParser* parser, int32_t offset, int32_t wherefrom)
{
    if (parser->buffer.byte_pos != -1) {
        /* The amount of data that has been read into the buffer is read ahead, so return it */
        offset -= (WAVBITBUFFER_BUFFER_SIZE - (parser->buffer.byte_pos + 1));
    }
    /* Move */
    fseek(parser->fp, offset, wherefrom);
    /* Clear the buffer */
    parser->buffer.byte_pos = -1;

    return WAV_ERROR_OK;
}

/* Discard the WAV file handle */
void WAV_Destroy(struct WAVFile* wavfile)
{
    uint32_t ch;

    /* Check for NULL and release */
#define NULLCHECK_AND_FREE(ptr) { \
    if ((ptr) != NULL) {            \
        free(ptr);                    \
        ptr = NULL;                   \
    }                               \
}

    if (wavfile != NULL) {
        for (ch = 0; ch < wavfile->format.num_channels; ch++) {
            NULLCHECK_AND_FREE(wavfile->data[ch]);
        }
        NULLCHECK_AND_FREE(wavfile->data);
        free(wavfile);
    }

#undef NULLCHECK_AND_FREE
}

/* Use a writer to output the header according to the file format */
static WAVError WAVWriter_PutWAVHeader(
        struct WAVWriter* writer, const struct WAVFileFormat* format)
{
    uint32_t filesize, pcm_data_size;

    /* Argument check */
    if (writer == NULL || format == NULL) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /* Format check */
    /* Only PCM is supported */
    if (format->data_format != WAV_DATA_FORMAT_PCM) {
        return WAV_ERROR_INVALID_FORMAT;
    }

    /* PCM data size */
    pcm_data_size
        = format->num_samples * (format->bits_per_sample / 8) * format->num_channels;

    /* File size */
    /* 44 is the number of bytes in the field from "RIFF" to (the size of "data") (not including any extensions) */
    filesize = pcm_data_size + 44;

    /* Output headers 'R', 'I', 'F', 'F' */
    if (WAVWriter_PutBits(writer, 'R', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'I', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'F', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'F', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* File size - 8 (size after this element) */
    if (WAVWriter_PutLittleEndianBytes(writer, 4, filesize - 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Output headers 'W', 'A', 'V', 'E' */
    if (WAVWriter_PutBits(writer, 'W', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'A', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'V', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'E', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Output fmt chunk headers 'f', 'm', 't', ' ' */
    if (WAVWriter_PutBits(writer, 'f', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'm', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 't', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, ' ', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Output the number of bytes in the fmt chunk (Note) Currently, it is set to 16 bytes */
    if (WAVWriter_PutLittleEndianBytes(writer, 4, 16) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Output format ID (Note) Currently set to 1 (Linear PCM) */
    if (WAVWriter_PutLittleEndianBytes(writer, 2,  1) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Number of channels */
    if (WAVWriter_PutLittleEndianBytes(writer, 2, format->num_channels) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Sampling rate */
    if (WAVWriter_PutLittleEndianBytes(writer, 4, format->sampling_rate) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Data speed (bytes/sec) */
    if (WAVWriter_PutLittleEndianBytes(writer, 4,
                format->sampling_rate * (format->bits_per_sample / 8) * format->num_channels)
            != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Size number per block */
    if (WAVWriter_PutLittleEndianBytes(writer, 2,
                (format->bits_per_sample / 8) * format->num_channels)
            != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    /* Quantization bit depth (bits per sample) */
    if (WAVWriter_PutLittleEndianBytes(writer, 2, format->bits_per_sample) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Output header of "data" chunk */
    if (WAVWriter_PutBits(writer, 'd', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'a', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 't', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };
    if (WAVWriter_PutBits(writer, 'a', 8) != WAV_ERROR_OK) { return WAV_ERROR_IO; };

    /* Number of bytes of waveform data */
    if (WAVWriter_PutLittleEndianBytes(writer, 4, pcm_data_size) != WAV_ERROR_OK) { return WAV_ERROR_IO; }

    return WAV_ERROR_OK;
}

/*
/* Write in little endian format
* Note) data may be swapped */
*/
static size_t WAVWrite_FWriteLittleEndian(
        void *data, size_t size, size_t ndata, FILE *fp)
{
    int x = 1;
    uint8_t *buffer;
    uint32_t i;

    /* In little endian environment, just fwrite */
    if ((size == 1) || (*((char *)&x) == 1)) {
        return fwrite(data, size, ndata, fp);
    }

    /* In big endian environment, rearrange before writing */
    buffer = (uint8_t *)data;

    switch (size) {
    case 2:
        for (i = 0; i < ndata; i++) {
            uint8_t a = buffer[2 * i];
            buffer[2 * i + 0] = buffer[2 * i + 1];
            buffer[2 * i + 1] = a;
        }
        break;
    case 3:
        for (i = 0; i < ndata; i++) {
            uint8_t a = buffer[3 * i];
            buffer[3 * i + 0] = buffer[3 * i + 2];
            buffer[3 * i + 2] = a;
        }
        break;
    case 4:
        for (i = 0; i < ndata; i++) {
            uint8_t a = buffer[4 * i];
            uint8_t b = buffer[4 * i + 1];
            buffer[4 * i + 0] = buffer[4 * i + 3];
            buffer[4 * i + 1] = buffer[4 * i + 2];
            buffer[4 * i + 2] = b;
            buffer[4 * i + 3] = a;
        }
        break;
    default:
        return 0;
    }

    return fwrite(data, size, ndata, fp);
}

/* Output PCM data using the writer */
static WAVError WAVWriter_PutWAVPcmData(
        struct WAVWriter* writer, const struct WAVFile* wavfile)
{
    uint32_t ch, smpl, progress;

    /* The buffer is empty */
    WAVWriter_Flush(writer);

    /* Write with channel interleaving */
    switch (wavfile->format.bits_per_sample) {
    case 8:
        {
            uint8_t *buffer;
            const uint32_t num_output_smpls_per_buffer = WAVBITBUFFER_BUFFER_SIZE / (sizeof(uint8_t) * wavfile->format.num_channels);
            progress = 0;
            while (progress < wavfile->format.num_samples) {
                const uint32_t num_process_smpls = WAV_Min(num_output_smpls_per_buffer, wavfile->format.num_samples - progress);
                const uint32_t num_output_smpls = num_process_smpls * wavfile->format.num_channels;
                buffer = (uint8_t *)writer->buffer.bytes;
                for (smpl = 0; smpl < num_process_smpls; smpl++) {
                    for (ch = 0; ch < wavfile->format.num_channels; ch++) {
                        (*buffer++) = (uint8_t)((WAVFile_PCM(wavfile, progress + smpl, ch) + 128) & 0xFF);
                    }
                }
                if (WAVWrite_FWriteLittleEndian(writer->buffer.bytes,
                            sizeof(uint8_t), num_output_smpls, writer->fp) < num_output_smpls) {
                    return WAV_ERROR_IO;
                }
                progress += num_process_smpls;
            }
        }
        break;
    case 16:
        {
            int16_t *buffer;
            const uint32_t num_output_smpls_per_buffer = (uint32_t)(WAVBITBUFFER_BUFFER_SIZE / (sizeof(int16_t) * wavfile->format.num_channels));
            progress = 0;
            while (progress < wavfile->format.num_samples) {
                const uint32_t num_process_smpls = WAV_Min(num_output_smpls_per_buffer, wavfile->format.num_samples - progress);
                const uint32_t num_output_smpls = num_process_smpls * wavfile->format.num_channels;
                buffer = (int16_t *)writer->buffer.bytes;
                for (smpl = 0; smpl < num_process_smpls; smpl++) {
                    for (ch = 0; ch < wavfile->format.num_channels; ch++) {
                        (*buffer++) = (int16_t)(WAVFile_PCM(wavfile, progress + smpl, ch) & 0xFFFF);
                    }
                }
                if (WAVWrite_FWriteLittleEndian(writer->buffer.bytes,
                            sizeof(int16_t), num_output_smpls, writer->fp) < num_output_smpls) {
                    return WAV_ERROR_IO;
                }
                progress += num_process_smpls;
            }
        }
        break;
    case 24:
        {
            uint8_t *buffer;
            const size_t int24_size = 3 * sizeof(uint8_t);
            const uint32_t num_output_smpls_per_buffer = (uint32_t)(WAVBITBUFFER_BUFFER_SIZE / (int24_size * wavfile->format.num_channels));
            progress = 0;
            while (progress < wavfile->format.num_samples) {
                const uint32_t num_process_smpls = WAV_Min(num_output_smpls_per_buffer, wavfile->format.num_samples - progress);
                const uint32_t num_output_smpls = num_process_smpls * wavfile->format.num_channels;
                const size_t output_size = num_output_smpls * int24_size;
                buffer = (uint8_t *)writer->buffer.bytes;
                for (smpl = 0; smpl < num_process_smpls; smpl++) {
                    for (ch = 0; ch < wavfile->format.num_channels; ch++) {
                        int32_t pcm = WAVFile_PCM(wavfile, progress + smpl, ch);
                        (*buffer++) = (uint8_t)((pcm >>  0) & 0xFF);
                        (*buffer++) = (uint8_t)((pcm >>  8) & 0xFF);
                        (*buffer++) = (uint8_t)((pcm >> 16) & 0xFF);
                    }
                }
                if (WAVWrite_FWriteLittleEndian(writer->buffer.bytes,
                            sizeof(uint8_t), output_size, writer->fp) < output_size) {
                    return WAV_ERROR_IO;
                }
                progress += num_process_smpls;
            }
        }
        break;
    case 32:
        {
            int32_t *buffer;
            const uint32_t num_output_smpls_per_buffer = (uint32_t)(WAVBITBUFFER_BUFFER_SIZE / (sizeof(int32_t) * wavfile->format.num_channels));
            progress = 0;
            while (progress < wavfile->format.num_samples) {
                const uint32_t num_process_smpls = WAV_Min(num_output_smpls_per_buffer, wavfile->format.num_samples - progress);
                const uint32_t num_output_smpls = num_process_smpls * wavfile->format.num_channels;
                buffer = (int32_t *)writer->buffer.bytes;
                for (smpl = 0; smpl < num_process_smpls; smpl++) {
                    for (ch = 0; ch < wavfile->format.num_channels; ch++) {
                        (*buffer++) = WAVFile_PCM(wavfile, progress + smpl, ch);
                    }
                }
                if (WAVWrite_FWriteLittleEndian(writer->buffer.bytes,
                            sizeof(int32_t), num_output_smpls, writer->fp) < num_output_smpls) {
                    return WAV_ERROR_IO;
                }
                progress += num_process_smpls;
            }
        }
        break;
    default:
        /* fprintf(stderr, "Unsupported bits per smpl format(=%d). \n", wavfile->format.bits_per_smpl); */
        return WAV_ERROR_INVALID_FORMAT;
    }


    return WAV_ERROR_OK;
}

/* Write file */
WAVApiResult WAV_WriteToFile(
        const char* filename, const struct WAVFile* wavfile)
{
    struct WAVWriter  writer;
    FILE*             fp;

    /* Argument check */
    if (filename == NULL || wavfile == NULL) {
        return WAV_APIRESULT_INVALID_PARAMETER;
    }

    /* Open wav file */
    fp = fopen(filename, "wb");
    if (fp == NULL) {
        /* fprintf(stderr, "Failed to open %s. \n", filename); */
        return WAV_APIRESULT_NG;
    }

    /* Writer initialization */
    WAVWriter_Initialize(&writer, fp);

    /* Write header */
    if (WAVWriter_PutWAVHeader(&writer, &wavfile->format) != WAV_ERROR_OK) {
        return WAV_APIRESULT_NG;
    }

    /* Write data */
    if (WAVWriter_PutWAVPcmData(&writer, wavfile) != WAV_ERROR_OK) {
        return WAV_APIRESULT_NG;
    }

    /* Writer finished */
    WAVWriter_Finalize(&writer);

    /* Close the file */
    fclose(fp);

    /* normal termination */
    return WAV_APIRESULT_OK;
}

/* Writer initialization */
static void WAVWriter_Initialize(struct WAVWriter* writer, FILE* fp)
{
    writer->fp                = fp;
    writer->bit_count         = 8;
    writer->bit_buffer        = 0;
    memset(&writer->buffer, 0, sizeof(struct WAVBitBuffer));
    writer->buffer.byte_pos   = 0;
}

/* Writer termination */
static void WAVWriter_Finalize(struct WAVWriter* writer)
{
    /* Write out any data remaining in the buffer */
    WAVWriter_Flush(writer);

    /* Clear members */
    writer->fp              = NULL;
    writer->bit_count       = 8;
    writer->bit_buffer      = 0;
    memset(&writer->buffer, 0, sizeof(struct WAVBitBuffer));
    writer->buffer.byte_pos = 0;
}

/* Write the lowest n_bits of val (in big endian) */
static WAVError WAVWriter_PutBits(struct WAVWriter* writer, uint64_t val, uint32_t n_bits)
{
    /* Invalid argument */
    if (writer == NULL) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /*
/* Sequentially output the most significant bits of val.
* In the first loop, fill in the remainder (number of bits required for output) before outputting.
* From the second loop onwards, output in 8-bit units. */
*/
    while (n_bits >= writer->bit_count) {
        n_bits -= writer->bit_count;
        writer->bit_buffer |= (uint8_t)WAV_GetLowerBits(writer->bit_count, val >> n_bits);

        /* Append to the buffer */
        writer->buffer.bytes[writer->buffer.byte_pos++] = (uint8_t)(writer->bit_buffer & 0xFF);

        /* Write out when buffer is full */
        if (writer->buffer.byte_pos == WAVBITBUFFER_BUFFER_SIZE) {
            if (fwrite(writer->buffer.bytes,
                        sizeof(uint8_t), WAVBITBUFFER_BUFFER_SIZE,
                        writer->fp) < WAVBITBUFFER_BUFFER_SIZE) {
                return WAV_ERROR_IO;
            }
            /* Reset write position */
            writer->buffer.byte_pos = 0;
        }

        writer->bit_buffer  = 0;
        writer->bit_count   = 8;
    }

    /*
/* Processing of fractional bits:
* Set the remaining bit to the upper bit of the buffer */
*/
    writer->bit_count -= n_bits;
    writer->bit_buffer |= (uint8_t)(WAV_GetLowerBits(n_bits, (uint32_t)val) << writer->bit_count);

    return WAV_ERROR_OK;
}

/* Output bit pattern in little endian */
static WAVError WAVWriter_PutLittleEndianBytes(
        struct WAVWriter* writer, uint32_t nbytes, uint64_t data)
{
    uint64_t out;
    uint32_t i_byte;

    /* Sort into little endian */
    out = 0;
    for (i_byte = 0; i_byte < nbytes; i_byte++) {
        out |= ((data >> (8 * (nbytes - i_byte - 1))) & 0xFFUL) << (8 * i_byte);
    }

    /* output */
    if (WAVWriter_PutBits(writer, out, (uint8_t)(nbytes * 8)) != WAV_ERROR_OK) {
        return WAV_ERROR_IO;
    }

    return WAV_ERROR_OK;
}

/* Clear the bits accumulated in the buffer */
static WAVError WAVWriter_Flush(struct WAVWriter* writer)
{
    /* Argument check */
    if (writer == NULL) {
        return WAV_ERROR_INVALID_PARAMETER;
    }

    /* Force output of remaining bits */
    if (writer->bit_count != 8) {
        if (WAVWriter_PutBits(writer, 0, (uint8_t)writer->bit_count) != WAV_ERROR_OK) {
            return WAV_ERROR_IO;
        }
        writer->bit_buffer = 0;
        writer->bit_count  = 8;
    }

    /* Flush any data remaining in the buffer */
    if (fwrite(writer->buffer.bytes,
                sizeof(uint8_t), (uint32_t)writer->buffer.byte_pos,
                writer->fp) < (size_t)writer->buffer.byte_pos) {
        return WAV_ERROR_IO;
    }
    /* Buffer remaining capacity is set to 0 */
    writer->buffer.byte_pos = 0;

    return WAV_ERROR_OK;
}

/* Get bit pattern in little endian */
static WAVError WAVParser_GetLittleEndianBytes(
        struct WAVParser* parser, uint32_t nbytes, uint64_t* bitsbuf)
{
    uint64_t tmp, ret;
    uint32_t i_byte;

    /* Get in big endian */
    if (WAVParser_GetBits(parser, nbytes * 8, &tmp) != WAV_ERROR_OK) {
        return WAV_ERROR_IO;
    }

    /* Sort into little endian */
    ret = 0;
    for (i_byte = 0; i_byte < nbytes; i_byte++) {
        ret |= ((tmp >> (8 * (nbytes - i_byte - 1))) & 0xFFUL) << (8 * i_byte);
    }
    *bitsbuf = ret;

    return WAV_ERROR_OK;
}

/* Get string using parser */
static WAVError WAVParser_GetString(
        struct WAVParser* parser, char* string_buffer, uint32_t string_length)
{
    uint32_t i_byte;
    uint64_t bitsbuf;

    assert(parser != NULL && string_buffer != NULL);

    /* Get string */
    for (i_byte = 0; i_byte < string_length; i_byte++) {
        /* Get 1 character */
        if (WAVParser_GetBits(parser, 8, &bitsbuf) != WAV_ERROR_OK) {
            return WAV_ERROR_IO;
        }
        string_buffer[i_byte] = (char)bitsbuf;
    }

    return WAV_ERROR_OK;
}

/* Use parser to get string/check for match */
static WAVError WAVParser_CheckSignatureString(
        struct WAVParser* parser, const char* signature, uint32_t signature_length)
{
    uint32_t i_byte;
    uint64_t bitsbuf;

    assert(parser != NULL && signature != NULL);

    /* Get/check string */
    for (i_byte = 0; i_byte < signature_length; i_byte++) {
        /* Get 1 character */
        if (WAVParser_GetBits(parser, 8, &bitsbuf) != WAV_ERROR_OK) {
            return WAV_ERROR_IO;
        }
        /* Signature check */
        if (signature[i_byte] != (char)bitsbuf) {
            /* fprintf(stderr, "Failed to check %s header signature. \n", signature); */
            return WAV_ERROR_INVALID_FORMAT;
        }
    }

    return WAV_ERROR_OK;
}
