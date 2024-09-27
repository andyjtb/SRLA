#include "srla_player.h"
#include <srla_decoder.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* Output request callback */
static void SRLAPlayer_SampleRequestCallback(int32_t **buffer, uint32_t num_channels, uint32_t num_samples);
/* Termination process */
static void exit_srla_player(void);

/* Global variables for playback control */
static struct SRLAHeader header = { 0, };
static uint32_t output_samples = 0;
static int32_t *decode_buffer[SRLA_MAX_NUM_CHANNELS] = { NULL, };
static uint32_t num_buffered_samples = 0;
static uint32_t buffer_pos = 0;
static uint32_t data_size = 0;
static uint8_t *data = NULL;
static uint32_t decode_offset = 0;
static struct SRLADecoder* decoder = NULL;

/* Main entry */
int main(int argc, char **argv)
{
    uint32_t i;
    SRLAApiResult ret;
    struct SRLADecoderConfig decoder_config;
    struct SRLAPlayerConfig player_config;

    /* Check arguments. If incorrect, show usage instructions. */
    if (argc != 2) {
        printf("Usage: %s SRLAFILE \n", argv[0]);
        return 1;
    }

    /* Load lnn file */
    {
        struct stat fstat;
        FILE* fp;
        const char *filename = argv[1];

        /* Open file */
        if ((fp = fopen(filename, "rb")) == NULL) {
            fprintf(stderr, "Failed to open %s \n", filename);
            return 1;
        }

        /* Get input file size / Allocate buffer space */
        stat(filename, &fstat);
        data_size = (uint32_t)fstat.st_size;
        data = (uint8_t *)malloc(data_size);

        /* Load data into the buffer area */
        if (fread(data, sizeof(uint8_t), data_size, fp) < data_size) {
            fprintf(stderr, "Failed to load %s data \n", filename);
            return 1;
        }

        fclose(fp);
    }

    /* Header decode */
    if ((ret = SRLADecoder_DecodeHeader(data, data_size, &header)) != SRLA_APIRESULT_OK) {
        fprintf(stderr, "Failed to get header information: %d \n", ret);
        return 1;
    }

    /* Create a decoder handle */
    decoder_config.max_num_channels = header.num_channels;
    decoder_config.max_num_parameters = SRLA_MAX_COEFFICIENT_ORDER;
    decoder_config.check_checksum = 1;
    if ((decoder = SRLADecoder_Create(&decoder_config, NULL, 0)) == NULL) {
        fprintf(stderr, "Failed to create decoder handle. \n");
        return 1;
    }

    /* Set the header in the decoder */
    if ((ret = SRLADecoder_SetHeader(decoder, &header)) != SRLA_APIRESULT_OK) {
        fprintf(stderr, "Failed to set header to decoder. \n");
        return 1;
    }

    /* Decode output area allocation */
    for (i = 0; i < header.num_channels; i++) {
        decode_buffer[i] = (int32_t *)malloc(sizeof(int32_t) * header.max_num_samples_per_block);
        memset(decode_buffer[i], 0, sizeof(int32_t) * header.max_num_samples_per_block);
    }

    /* Advance the decoding position by the header */
    decode_offset = SRLA_HEADER_SIZE;

    /* Player initialization */
    player_config.sampling_rate = header.sampling_rate;
    player_config.num_channels = header.num_channels;
    player_config.bits_per_sample = header.bits_per_sample;
    player_config.sample_request_callback = SRLAPlayer_SampleRequestCallback;
    SRLAPlayer_Initialize(&player_config);

    /* After this, proceed according to the callback request */
    while (1) { ; }

    return 0;
}

/* Output request callback */
static void SRLAPlayer_SampleRequestCallback(int32_t **buffer, uint32_t num_channels, uint32_t num_samples)
{
    uint32_t ch, smpl;

    for (smpl = 0; smpl < num_samples; smpl++) {
        /* Decode immediately when buffer is full */
        if (buffer_pos >= num_buffered_samples) {
            uint32_t decode_size;
            if (SRLADecoder_DecodeBlock(decoder,
                        &data[decode_offset], data_size - decode_offset,
                        decode_buffer, header.num_channels, header.max_num_samples_per_block,
                        &decode_size, &num_buffered_samples) != SRLA_APIRESULT_OK) {
                fprintf(stderr, "decoding error! \n");
                exit(1);
            }
            buffer_pos = 0;
            decode_offset += decode_size;
        }

        /* Copy to output buffer area */
        for (ch = 0; ch < num_channels; ch++) {
            buffer[ch][smpl] = decode_buffer[ch][buffer_pos];
        }
        buffer_pos++;
        output_samples++;

        /* End processing as soon as playback ends */
        if (output_samples >= header.num_samples) {
            exit_srla_player();
        }
    }

    /* Progress display */
    printf("playing... %7.3f / %7.3f \r",
            (double)output_samples / header.sampling_rate, (double)header.num_samples / header.sampling_rate);
    fflush(stdout);
}

/* Termination process */
static void exit_srla_player(void)
{
    uint32_t i;

    SRLAPlayer_Finalize();

    for (i = 0; i < header.num_channels; i++) {
        free(decode_buffer[i]);
    }
    SRLADecoder_Destroy(decoder);
    free(data);

    exit(0);
}
