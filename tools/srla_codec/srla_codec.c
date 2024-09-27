#include <srla_encoder.h>
#include <srla_decoder.h>
#include "wav.h"
#include "command_line_parser.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* Default preset */
#define DEFALUT_PRESET_INDEX 2
/* Default maximum block samples */
#define DEFALUT_MAX_NUM_BLOCK_SAMPLES 4096
/* Default variable block division number */
#define DEFALUT_NUM_VARIABLE_BLOCK_DIVISIONS 1
/* Maximum index of parameter preset */
#define SRLA_MAX_PARAMETER_PRESETS_INDEX 4
#if SRLA_MAX_PARAMETER_PRESETS_INDEX != (SRLA_NUM_PARAMETER_PRESETS - 1)
#error "Max parameter presets mismatched to number of parameter presets!"
#endif
/* Convert the macro content to a string */
#define PRE_TOSTRING(arg) #arg
#define TOSTRING(arg) PRE_TOSTRING(arg)

/* Select the smaller of a and b */
#define SRLACODEC_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* Command line specification */
static struct CommandLineParserSpecification command_line_spec[] = {
    { 'e', "encode", "Encode mode",
        COMMAND_LINE_PARSER_FALSE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'd', "decode", "Decode mode",
        COMMAND_LINE_PARSER_FALSE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'm', "mode", "Specify compress mode: 0(fast), ..., " TOSTRING(SRLA_MAX_PARAMETER_PRESETS_INDEX) "(high compression) (default:" TOSTRING(DEFALUT_PRESET_INDEX) ")",
        COMMAND_LINE_PARSER_TRUE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'B', "max-block-size", "Specify max number of block samples (default:" TOSTRING(DEFALUT_MAX_NUM_BLOCK_SAMPLES) ")",
        COMMAND_LINE_PARSER_TRUE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'V', "variable-block-divisions", "Specify number of variable block-size divisions (default:" TOSTRING(DEFALUT_NUM_VARIABLE_BLOCK_DIVISIONS) ")",
        COMMAND_LINE_PARSER_TRUE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'c', "no-checksum-check", "Whether to NOT check checksum at decoding (default:no)",
        COMMAND_LINE_PARSER_FALSE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'h', "help", "Show command help message",
        COMMAND_LINE_PARSER_FALSE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 'v', "version", "Show version information",
        COMMAND_LINE_PARSER_FALSE, NULL, COMMAND_LINE_PARSER_FALSE },
    { 0, }
};

/* Encode returns 0 on success, non-zero on failure */
static int do_encode(const char *in_filename, const char *out_filename,
    uint32_t encode_preset_no, uint32_t max_num_block_samples, uint32_t variable_block_num_divisions)
{
    FILE *out_fp;
    struct WAVFile *in_wav;
    struct SRLAEncoder *encoder;
    struct SRLAEncoderConfig config;
    struct SRLAEncodeParameter parameter;
    struct stat fstat;
    uint8_t *buffer;
    uint32_t buffer_size, encoded_data_size;
    SRLAApiResult ret;
    uint32_t ch, smpl, num_channels, num_samples;
    SRLAApiResult (*encode_function)(struct SRLAEncoder* encoder,
        const int32_t* const* input, uint32_t num_samples,
        uint8_t * data, uint32_t data_size, uint32_t * output_size);

    /* Create the encoder */
    config.max_num_channels = SRLA_MAX_NUM_CHANNELS;
    config.min_num_samples_per_block = max_num_block_samples >> variable_block_num_divisions;
    config.max_num_samples_per_block = max_num_block_samples;
    config.max_num_parameters = SRLA_MAX_COEFFICIENT_ORDER;
    if ((encoder = SRLAEncoder_Create(&config, NULL, 0)) == NULL) {
        fprintf(stderr, "Failed to create encoder handle. \n");
        return 1;
    }

    /* Open WAV file */
    if ((in_wav = WAV_CreateFromFile(in_filename)) == NULL) {
        fprintf(stderr, "Failed to open %s. \n", in_filename);
        return 1;
    }
    num_channels = in_wav->format.num_channels;
    num_samples = in_wav->format.num_samples;

    /* Encoding parameter set */
    parameter.num_channels = (uint16_t)num_channels;
    parameter.bits_per_sample = (uint16_t)in_wav->format.bits_per_sample;
    parameter.sampling_rate = in_wav->format.sampling_rate;
    parameter.min_num_samples_per_block = max_num_block_samples >> variable_block_num_divisions;
    parameter.max_num_samples_per_block = max_num_block_samples;
    /* Reflecting presets */
    parameter.preset = (uint8_t)encode_preset_no;
    if ((ret = SRLAEncoder_SetEncodeParameter(encoder, &parameter)) != SRLA_APIRESULT_OK) {
        fprintf(stderr, "Failed to set encode parameter: %d \n", ret);
        return 1;
    }

    /* Keep track of the size of the input file */
    stat(in_filename, &fstat);
    /* Assuming it will not be larger than twice the input wav size */
    buffer_size = (uint32_t)(2 * fstat.st_size);

    /* Create the encoded data area */
    buffer = (uint8_t *)malloc(buffer_size);

    /* Select the encoding function */
    encode_function = (variable_block_num_divisions == 0) ? SRLAEncoder_EncodeBlock : SRLAEncoder_EncodeOptimalPartitionedBlock;

    /* Execute encoding */
    {
        uint8_t *data_pos = buffer;
        uint32_t write_offset, progress;
        struct SRLAHeader header;

        write_offset = 0;

        /* Header encoding */
        header.num_channels = (uint16_t)num_channels;
        header.num_samples = num_samples;
        header.sampling_rate = parameter.sampling_rate;
        header.bits_per_sample = parameter.bits_per_sample;
        header.preset = parameter.preset;
        header.max_num_samples_per_block = parameter.max_num_samples_per_block;
        if ((ret = SRLAEncoder_EncodeHeader(&header, data_pos, buffer_size))
                != SRLA_APIRESULT_OK) {
            fprintf(stderr, "Failed to encode header! ret:%d \n", ret);
            return 1;
        }
        data_pos += SRLA_HEADER_SIZE;
        write_offset += SRLA_HEADER_SIZE;

        /* Encode blocks in chronological order */
        progress = 0;
        while (progress < num_samples) {
            uint32_t ch, write_size;
            const int32_t *input_ptr[SRLA_MAX_NUM_CHANNELS];
            /* Determine the number of samples to encode */
            const uint32_t num_encode_samples = SRLACODEC_MIN(max_num_block_samples, num_samples - progress);

            /* Set sample reference position */
            for (ch = 0; ch < (uint32_t)num_channels; ch++) {
                input_ptr[ch] = &(WAVFile_PCM(in_wav, progress, ch));
            }

            /* Block encoding */
            if ((ret = encode_function(encoder,
                    input_ptr, num_encode_samples,
                    data_pos, buffer_size - write_offset, &write_size)) != SRLA_APIRESULT_OK) {
                fprintf(stderr, "Failed to encode! ret:%d \n", ret);
                return 1;
            }

            /* Progress update */
            data_pos += write_size;
            write_offset += write_size;
            progress += num_encode_samples;

            /* Progress display */
            printf("progress... %5.2f%% \r", (double)((progress * 100.0) / num_samples));
            fflush(stdout);
        }

        /* Get write size */
        encoded_data_size = write_offset;
    }

    /* Write file */
    out_fp = fopen(out_filename, "wb");
    if (fwrite(buffer, sizeof(uint8_t), encoded_data_size, out_fp) < encoded_data_size) {
        fprintf(stderr, "File output error! %d \n", ret);
        return 1;
    }

    /* Display compression result summary */
    printf("finished: %d -> %d (%6.2f %%) \n",
            (uint32_t)fstat.st_size, encoded_data_size, (double)((100.0 * encoded_data_size) / (double)fstat.st_size));

    /* Dispose of resources */
    fclose(out_fp);
    free(buffer);
    WAV_Destroy(in_wav);
    SRLAEncoder_Destroy(encoder);

    return 0;
}

/* Returns 0 on success, non-zero on failure */
static int do_decode(const char *in_filename, const char *out_filename, uint8_t check_checksum)
{
    FILE* in_fp;
    struct WAVFile* out_wav;
    struct WAVFileFormat wav_format;
    struct stat fstat;
    struct SRLADecoder* decoder;
    struct SRLADecoderConfig config;
    struct SRLAHeader header;
    uint8_t* buffer;
    uint32_t ch, smpl, buffer_size;
    SRLAApiResult ret;

    /* Create a decoder handle */
    config.max_num_channels = SRLA_MAX_NUM_CHANNELS;
    config.max_num_parameters = SRLA_MAX_COEFFICIENT_ORDER;
    config.check_checksum = check_checksum;
    if ((decoder = SRLADecoder_Create(&config, NULL, 0)) == NULL) {
        fprintf(stderr, "Failed to create decoder handle. \n");
        return 1;
    }

    /* Open input file */
    in_fp = fopen(in_filename, "rb");
    /* Get input file size / Allocate buffer space */
    stat(in_filename, &fstat);
    buffer_size = (uint32_t)fstat.st_size;
    buffer = (uint8_t *)malloc(buffer_size);
    /* Load data into the buffer area */
    fread(buffer, sizeof(uint8_t), buffer_size, in_fp);
    fclose(in_fp);

    /* Header decode */
    if ((ret = SRLADecoder_DecodeHeader(buffer, buffer_size, &header))
            != SRLA_APIRESULT_OK) {
        fprintf(stderr, "Failed to get header information: %d \n", ret);
        return 1;
    }

    /* Generate output wav handle */
    wav_format.data_format     = WAV_DATA_FORMAT_PCM;
    wav_format.num_channels    = header.num_channels;
    wav_format.sampling_rate   = header.sampling_rate;
    wav_format.bits_per_sample = header.bits_per_sample;
    wav_format.num_samples     = header.num_samples;
    if ((out_wav = WAV_Create(&wav_format)) == NULL) {
        fprintf(stderr, "Failed to create wav handle. \n");
        return 1;
    }

    /* Bulk decode */
    if ((ret = SRLADecoder_DecodeWhole(decoder,
                    buffer, buffer_size,
                    (int32_t **)out_wav->data, out_wav->format.num_channels, out_wav->format.num_samples))
                != SRLA_APIRESULT_OK) {
        fprintf(stderr, "Decoding error! %d \n", ret);
        return 1;
    }

    /* Export WAV file */
    if (WAV_WriteToFile(out_filename, out_wav) != WAV_APIRESULT_OK) {
        fprintf(stderr, "Failed to write wav file. \n");
        return 1;
    }

    free(buffer);
    WAV_Destroy(out_wav);
    SRLADecoder_Destroy(decoder);

    return 0;
}

/* Show usage */
static void print_usage(char** argv)
{
    printf("Usage: %s [options] INPUT_FILE_NAME OUTPUT_FILE_NAME \n", argv[0]);
}

/* Display version information */
static void print_version_info(void)
{
    printf("SRLA -- SVR-FIR Lossless Audio codec Version.%d \n", SRLA_CODEC_VERSION);
}

/* Main entry */
int main(int argc, char** argv)
{
    const char *filename_ptr[2] = { NULL, NULL };
    const char *input_file;
    const char *output_file;

    /* Not enough arguments */
    if (argc == 1) {
        print_usage(argv);
        /* Prompts users to display help so that first-time users don't get stuck */
        printf("Type `%s -h` to display command helps. \n", argv[0]);
        return 1;
    }

    /* Command line parsing */
    if (CommandLineParser_ParseArguments(command_line_spec,
                argc, (const char *const *)argv, filename_ptr, sizeof(filename_ptr) / sizeof(filename_ptr[0]))
            != COMMAND_LINE_PARSER_RESULT_OK) {
        return 1;
    }

    /* Determine whether help or version information is displayed */
    if (CommandLineParser_GetOptionAcquired(command_line_spec, "help") == COMMAND_LINE_PARSER_TRUE) {
        print_usage(argv);
        printf("options: \n");
        CommandLineParser_PrintDescription(command_line_spec);
        return 0;
    } else if (CommandLineParser_GetOptionAcquired(command_line_spec, "version") == COMMAND_LINE_PARSER_TRUE) {
        print_version_info();
        return 0;
    }

    /* Get input file name */
    if ((input_file = filename_ptr[0]) == NULL) {
        fprintf(stderr, "%s: input file must be specified. \n", argv[0]);
        return 1;
    }

    /* Get output file name */
    if ((output_file = filename_ptr[1]) == NULL) {
        fprintf(stderr, "%s: output file must be specified. \n", argv[0]);
        return 1;
    }

    /* Encode and decode cannot be specified at the same time */
    if ((CommandLineParser_GetOptionAcquired(command_line_spec, "decode") == COMMAND_LINE_PARSER_TRUE)
            && (CommandLineParser_GetOptionAcquired(command_line_spec, "encode") == COMMAND_LINE_PARSER_TRUE)) {
        fprintf(stderr, "%s: encode and decode mode cannot specify simultaneously. \n", argv[0]);
        return 1;
    }

    if (CommandLineParser_GetOptionAcquired(command_line_spec, "decode") == COMMAND_LINE_PARSER_TRUE) {
        /* Decode */
        uint8_t crc_check = 1;
        /* Get CRC invalid flag */
        if (CommandLineParser_GetOptionAcquired(command_line_spec, "no-crc-check") == COMMAND_LINE_PARSER_TRUE) {
            crc_check = 0;
        }
        /* Execute batch decoding */
        if (do_decode(input_file, output_file, crc_check) != 0) {
            fprintf(stderr, "%s: failed to decode %s. \n", argv[0], input_file);
            return 1;
        }
    } else if (CommandLineParser_GetOptionAcquired(command_line_spec, "encode") == COMMAND_LINE_PARSER_TRUE) {
        /* Encode */
        uint32_t encode_preset_no = DEFALUT_PRESET_INDEX;
        uint32_t max_num_block_samples = DEFALUT_MAX_NUM_BLOCK_SAMPLES;
        uint32_t variable_block_num_divisions = DEFALUT_NUM_VARIABLE_BLOCK_DIVISIONS;
        /* Get encoding preset number */
        if (CommandLineParser_GetOptionAcquired(command_line_spec, "mode") == COMMAND_LINE_PARSER_TRUE) {
            char *e;
            const char *lstr = CommandLineParser_GetArgumentString(command_line_spec, "mode");
            encode_preset_no = (uint32_t)strtol(lstr, &e, 10);
            if (*e != '\0') {
                fprintf(stderr, "%s: invalid encode preset number. (irregular character found in %s at %s)\n", argv[0], lstr, e);
                return 1;
            }
            if (encode_preset_no >= SRLA_NUM_PARAMETER_PRESETS) {
                fprintf(stderr, "%s: encode preset number is out of range. \n", argv[0]);
                return 1;
            }
        }
        /* Get number of samples per block */
        if (CommandLineParser_GetOptionAcquired(command_line_spec, "max-block-size") == COMMAND_LINE_PARSER_TRUE) {
            char *e;
            const char *lstr = CommandLineParser_GetArgumentString(command_line_spec, "max-block-size");
            max_num_block_samples = (uint32_t)strtol(lstr, &e, 10);
            if (*e != '\0') {
                fprintf(stderr, "%s: invalid number of block samples. (irregular character found in %s at %s)\n", argv[0], lstr, e);
                return 1;
            }
            if ((max_num_block_samples == 0) || (max_num_block_samples >= (1U << 16))) {
                fprintf(stderr, "%s: number of block samples is out of range. \n", argv[0]);
                return 1;
            }
        }
        /* Variable block encoding division number */
        if (CommandLineParser_GetOptionAcquired(command_line_spec, "variable-block-divisions") == COMMAND_LINE_PARSER_TRUE) {
            char *e;
            const char *lstr = CommandLineParser_GetArgumentString(command_line_spec, "variable-block-divisions");
            variable_block_num_divisions = (uint32_t)strtol(lstr, &e, 10);
            if (*e != '\0') {
                fprintf(stderr, "%s: invalid number of variable block divisions. (irregular character found in %s at %s)\n", argv[0], lstr, e);
                return 1;
            }
            if ((max_num_block_samples >> variable_block_num_divisions) == 0) {
                fprintf(stderr, "%s: number of variable block divisions is too large. \n", argv[0]);
                return 1;
            }
        }
        /* Execute batch encoding */
        if (do_encode(input_file, output_file, encode_preset_no, max_num_block_samples, variable_block_num_divisions) != 0) {
            fprintf(stderr, "%s: failed to encode %s. \n", argv[0], input_file);
            return 1;
        }
    } else {
        fprintf(stderr, "%s: decode(-d) or encode(-e) option must be specified. \n", argv[0]);
        return 1;
    }

    return 0;
}
